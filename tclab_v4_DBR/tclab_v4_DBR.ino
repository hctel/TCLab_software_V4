/*
  TCLab Temperature Control Lab Firmware
  Jeffrey Kantor, Bill Tubbs, John Hedengren
  May 2019

  This firmware provides a high level interface to the Temperature Control Lab. The
  firmware scans the serial port for commands. Commands are case-insensitive. Any
  unrecognized command results in sleep model. Each command returns a result string.

  A         software restart. Returns "Start".
  LED float set LED to float for 10 sec. range 0 to 100. Returns actual float
  P1 float  set pwm limit on heater 1, range 0 to 255. Default 200. Returns P1.
  P2 float  set pwm limit on heater 2, range 0 to 255. Default 100. Returns P2.
  Q1 float  set Heater 1, range 0 to 100. Returns value of Q1.
  Q2 float  set Heater 2, range 0 to 100. Returns value of Q2.
  R1        get value of Heater 1, range 0 to 100
  R2        get value of Heater 2, range 0 to 100
  SCAN      get values T1 T2 Q1 Q1 in line delimited values
  T1        get Temperature T1. Returns value of T1 in °C.
  T2        get Temperature T2. Returns value of T2 in °C.
  VER       get firmware version string
  X         stop, enter sleep mode. Returns "Stop"

  Limits on the heater can be configured with the constants below.

  Status is indicated by LED1 on the Temperature Control Lab. Status conditions are:

      LED1        LED1
      Brightness  State
      ----------  -----
      dim         steady     Normal operation, heaters off
      bright      steady     Normal operation, heaters on
      dim         blinking   High temperature alarm on, heaters off
      bright      blinking   High temperature alarm on, heaters on

  The Temperature Control Lab shuts down the heaters if it receives no host commands
  during a timeout period (configure below), receives an "X" command, or receives
  an unrecognized command from the host.

  The constants can be used to configure the firmware.

  Version History
      1.0.1 first version included in the tclab package
      1.1.0 added R1 and R2 commands to read current heater values
            modified heater values to units of percent of full power
            added P1 and P2 commands to set heater power limits
            rewrote readCommand to avoid busy states
            simplified LED status model
      1.2.0 added LED command
      1.2.1 correctly reset heater values on close
            added version history
      1.2.2 shorten version string for better display by TCLab
      1.2.3 move baudrate to from 9600 to 115200
      1.3.0 add SCAN function
            report board type in version string
      1.4.0 changed Q1 and Q2 to float from int
      1.4.1 fix missing Serial.flush() at end of command loop
      1.4.2 fix bug with X command
      1.4.3 required Arduino IDE Version >= 1.0.0
      1.5.0 remove webusb
      1.6.0 average 10 temperature measurements to reduce noise
      1.7.0 change of various variables (DBR)
      1.8.0 change dynamics to give a S-shape (DBR)
      1.9.0 allow for load disturbances (DBR)
*/

#include "Arduino.h"

// determine board type
#if defined(__AVR_ATmega328P__) || defined(__AVR_ATmega168__)
  String boardType = "Arduino Uno";
#elif defined(__AVR_ATmega32U4__) || defined(__AVR_ATmega16U4__)
  String boardType = "Arduino Leonardo/Micro";
#elif defined(__AVR_ATmega1280__) || defined(__AVR_ATmega2560__)
  String boardType = "Arduino Mega";
#else 
  String boardType = "Unknown board";
#endif

// Enable debugging output
const bool DEBUG = false;

// Define Q1/Q2 filter constants
const float Tlag1 = 120; // First lag time constant (adjust as needed)
const float Tlag2 = 30; // Second lag time constant (adjust as needed)
const float Tlead = 150; // Lead time constant (adjust as needed)

// constants
const String vers = "1.9.0";   // version of this firmware
const long baud = 115200;      // serial baud rate
const char sp = ' ';           // command separator
const char nl = '\n';          // command terminator

// pin numbers corresponding to signals on the TC Lab Shield
const int pinT1   = 0;         // T1
const int pinT2   = 2;         // T2
const int pinQ1   = 3;         // Q1
const int pinQ2   = 5;         // Q2
const int pinLED1 = 9;         // LED1

// temperature alarm limits
const int limT1   = 50;       // T1 high alarm (°C)
const int limT2   = 50;       // T2 high alarm (°C)

// LED1 levels
const int hiLED   =  60;       // hi LED
const int loLED   = hiLED/16;  // lo LED

// global variables
unsigned long lastTimeQ1 = 0; // Variable to hold the last loop time
unsigned long lastTimeQ2 = 0; // Variable to hold the last loop time
char Buffer[64];               // buffer for parsing serial input
int buffer_index = 0;          // index for Buffer
String cmd;                    // command
float val;                     // command value
int ledStatus;                 // 1: loLED
                               // 2: hiLED
                               // 3: loLED blink
                               // 4: hiLED blink
long ledTimeout = 0;           // when to return LED to normal operation
float LED = 100;               // LED override brightness
float P1 = 220;                // heater 1 power limit in units of pwm. Range 0 to 255
float P2 = 220;                // heater 2 power limit in units in pwm, Range 0 to 255
float Q1 = 0;                  // last value written to heater 1 in units of percent
boolean DL = false;            // DL = true ==> Load disturbance  
float MulQ1 = 1;               // multiplier for Q1 used for load disturbance
float Q2 = 0;                  // last value written to heater 2 in units of percent
int alarmStatus;               // hi temperature alarm status
boolean newData = false;       // boolean flag indicating new command

int n = 10;                    // number of samples for each temperature measurement

float QfLL1 = 0.0; // Initial value for QfLL1
float QfLL2 = 0.0; // Initial value for QfLL2
float Qp1 = 0.0; // Previous value for Q1
float Qp2 = 0.0; // Previous value for Q2

void readCommand() {
  while (Serial && (Serial.available() > 0) && (newData == false)) {
    int byte = Serial.read();
    if ((byte != '\r') && (byte != nl) && (buffer_index < 64)) {
      Buffer[buffer_index] = byte;
      buffer_index++;
    }
    else {
      newData = true;
    }
  }   
}

// for debugging with the serial monitor in Arduino IDE
void echoCommand() {
  if (newData) {
    Serial.write("Received Command: ");
    Serial.write(Buffer, buffer_index);
    Serial.write(nl);
    Serial.flush();
  }
}

// return thermister temperature in °C
inline float readTemperature(int pin) {
  float degC = 0.;
  for (int i = 0; i < n; i++) {
    // (analogRead * 3300.0/1024.0 - 500.0)/10.0
    degC += analogRead(pin) * 0.322265625 - 50.0;
  }
  // return average from n reads
  return degC / n;
}

void parseCommand(void) {
  if (newData) {
    String read_ = String(Buffer);

    // separate command from associated data
    int idx = read_.indexOf(sp);
    cmd = read_.substring(0, idx);
    cmd.trim();
    cmd.toUpperCase();

    // extract data. toFloat() returns 0 on error
    String data = read_.substring(idx + 1);
    data.trim();
    val = data.toFloat();

    // reset parameter for next command
    memset(Buffer, 0, sizeof(Buffer));
    buffer_index = 0;
    newData = false;
  }
}

void sendResponse(String msg) {
  Serial.println(msg);
}

void dispatchCommand(void) {
  if (cmd == "A") {
    setHeater1(0);
    setHeater2(0);
    sendResponse("Start");
  }
  else if (cmd == "LED") {
    ledTimeout = millis() + 10000;
    LED = max(0, min(100, val));
    sendResponse(String(LED));
  }
  else if (cmd == "P1") {
    P1 = max(0, min(255, val));
    sendResponse(String(P1));
  }
  else if (cmd == "P2") {
    P2 = max(0, min(255, val));
    sendResponse(String(P2));
  }
  else if (cmd == "Q1") {
    setHeater1(val);    
    sendResponse(String(Q1));
  }
  else if (cmd == "Q2") {
    setHeater2(val);    
    sendResponse(String(Q2));
  }
  else if (cmd == "DL") {
    if (val == true) {
      MulQ1 = 0.7;
    } else {
      MulQ1 = 1.0;
    }
    sendResponse(String(DL));    
  }
  else if (cmd == "RMULQ1") {
    sendResponse(String(MulQ1));
  }  
  else if (cmd == "R1") {
    sendResponse(String(Q1));
  }
  else if (cmd == "R2") {
    sendResponse(String(Q2));
  }
  else if (cmd == "SCAN") {
    sendResponse(String(readTemperature(pinT1)));
    sendResponse(String(readTemperature(pinT2)));
    sendResponse(String(Q1));
    sendResponse(String(Q2));
  }
  else if (cmd == "T1") {
    sendResponse(String(readTemperature(pinT1)));
  }
  else if (cmd == "T2") {
    sendResponse(String(readTemperature(pinT2)));
  }
  else if (cmd == "VER") {
    sendResponse("TCLab Firmware " + vers + " " + boardType);
  }
  else if (cmd == "X") {
    setHeater1(0);
    setHeater2(0);
  }
  else if (cmd.length() > 0) {
    setHeater1(0);
    setHeater2(0);
    sendResponse(cmd);
  }
  Serial.flush();
  cmd = "";
}

void checkAlarm(void) {
  if ((readTemperature(pinT1) > limT1) or (readTemperature(pinT2) > limT2)) {
    alarmStatus = 1;
  }
  else {
    alarmStatus = 0;
  }
}

void updateStatus(void) {
  // determine led status
  ledStatus = 1;
  if ((Q1 > 0) or (Q2 > 0)) {
    ledStatus = 2;
  }
  if (alarmStatus > 0) {
    ledStatus += 2;
  }
  // update led depending on ledStatus
  if (millis() < ledTimeout) {        // override led operation
    analogWrite(pinLED1, LED);
  }
  else {
    switch (ledStatus) {
      case 1:  // normal operation, heaters off
        analogWrite(pinLED1, loLED);
        break;
      case 2:  // normal operation, heater on
        analogWrite(pinLED1, hiLED);
        break;
      case 3:  // high temperature alarm, heater off
        if ((millis() % 2000) > 1000) {
          analogWrite(pinLED1, loLED);
        } else {
          analogWrite(pinLED1, loLED/4);
        }
        break;
      case 4:  // high temperature alarm, heater on
        if ((millis() % 2000) > 1000) {
          analogWrite(pinLED1, hiLED);
        } else {
          analogWrite(pinLED1, loLED);
        }
        break;
    }   
  }
}

// set Heater 1
void setHeater1(float qval) {
  Q1 = max(0., min(qval, 100.));
}

// set Heater 2
void setHeater2(float qval) {
  Q2 = max(0., min(qval, 100.));
}

// Function to calculate sampling time
float calculateSamplingTime(unsigned long& lastTime) {
  unsigned long currentTime = millis();
  float Ts = (currentTime - lastTime) / 1000.0; // Convert to seconds
  lastTime = currentTime;
  return Ts;
}

// Function to perform filtering
float Filter(float Q, float& QfLL, float& Qp, unsigned long& lastTime) {
  
  // Calculate sampling time
  float Ts = calculateSamplingTime(lastTime);
  
  // Static variables initialized once
  static float Qfilt = 0.0; 
  static float KpLL = 0.0; 
  static float KpFO = 0.0; 
  
  KpLL = Ts/Tlag1;
  KpFO = Ts/Tlag2;
  
  // Perform filtering
  QfLL = (1/(1+KpLL)) * QfLL + (KpLL/(1+KpLL)) * ((1+Tlead/Ts)*Q - (Tlead/Ts)*Qp);
  Qfilt = (1/(1+KpFO)) * Qfilt + (KpFO/(1+KpFO)) * QfLL;
  Qp = Q;
  
  return Qfilt;
}

// Heater 1 filter
void Heater1Filter() {
  float filteredQ1 = Filter(Q1, QfLL1, Qp1, lastTimeQ1);
  filteredQ1 = max(0., min(filteredQ1, 100.));
  analogWrite(pinQ1, (MulQ1*filteredQ1*P1)/100);
}

// Heater 2 filter
void Heater2Filter() {
  float filteredQ2 = Filter(Q2, QfLL2, Qp2, lastTimeQ2);
  filteredQ2 = max(0., min(filteredQ2, 100.));
  analogWrite(pinQ2, (filteredQ2*P2)/100);
}

// arduino startup
void setup() {
  analogReference(EXTERNAL);
  while (!Serial) {
    ; // wait for serial port to connect.
  }
  Serial.begin(baud);
  Serial.flush();
  setHeater1(0);
  setHeater2(0);
  ledTimeout = millis() + 1000;
}

// arduino main event loop
void loop() {
  readCommand();
  if (DEBUG) echoCommand();
  parseCommand();
  dispatchCommand();
  Heater1Filter();
  Heater2Filter();
  checkAlarm();
  updateStatus();
}
