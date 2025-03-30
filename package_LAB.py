import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

CLP_FF = 0
CLP_NOFF = 1
OLP_FF = 2
OLP_NOFF = 3
CLP_FF_FOLLOW = 4

def LL_RT(MV,Kp,Tlead,Tlag,Ts,PV,PVInit=0,method='EBD'):
    
    """
    The function "LL_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :Tlead: lead time constant [s]
    :Tlag: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method (Not implemented yet)
    
    The function "LL_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    
    if (Tlag != 0):
        K = Ts/Tlag
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if len(MV) == 0:
                MV.append(0)
            elif len(MV) == 1:
                MV.append(0)
            else:
                if method == 'EBD':
                    PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*((1+Tlead/Tlag)*MV[-1] - Tlead/Tlag*MV[-2]))
                elif method == 'EFD':
                    PV.append((1-K)*PV[-1] + K*Kp*((1+Tlead/Tlag)*MV[-1] - Tlead/Tlag*MV[-2]))
                    #TODO: Implement TRAP method
                else:
                    PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*((1+Tlead/Tlag)*MV[-1] - Tlead/Tlag*MV[-2]))
    else:
        PV.append(Kp*MV[-1])

#-----------------------------------
def LL(MV,Kp,Tlead, Tlag,Ts,MVInit=0,PVInit=0,method='EBD'):
    
    """
    The function "FOPDT" DOES NOT need to be included in a "for or while loop": this block is for offline use.
    
    :MV: input vector
    :Kp: process gain
    :T: lag time constant [s]
    :theta: delay [s]
    :Ts: sampling period [s]
    :MVInit: (optional: default value is 0)    
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
        
    :return: simulated FOPDT output vector         
    
    The function "FOPDT" returns the simulated output FOPDT vector from the input vector "MV" and the input parameters.
    """    
    
    MVDelay = []
    MVTemp = []
    PVSim = []
    
    for i in range(0,len(MV)):
        MVTemp.append(MV[i])
        LL_RT(MVDelay,Kp,Tlead,Tlag,Ts,PVSim,PVInit,method)
            
    return PVSim

#-----------------------------------

def PID_RT(SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E, ManFF=False, PVInit=0, method='EBD-EBD'):

    """
    The function "PID_RT" needs to be included in a -for or while loop".

    :SP: SP (or SetPoint) vector
    :PV: PV (or Process Value) vector 
    :Man: Man (or Manual controller mode) vector (True or False) 
    :MVMan: MVMan (or Manual value for MV) vector 
    :MVFF: MVFF (or Feedforward) vector 

    :Kc: controller gain 
    :Ti: integral time constant [s] 
    :Td: derivative time constant [s] 
    :alpha: Tfd = alpha*Td where Tfd is the derivative filter time constant [s] 
    :Ts: sampling period [s]
    
    :MVMin: minimum value for MV (used for saturation and anti wind-up) 
    :MVMax: maximum value for MV (used for saturation and anti wind-up) 

    :MV: MV (or Manipulated Value) vector 
    :MVP: MVP (or Propotional part of MV) vector 
    :MVI: MVI (or Integral part of MV) vector 
    :MVD: MVD (or Derivative part of MV) vector 
    :E: E (or control Error) vector 

    :ManFF: Activated FF in manual mode (optional: default boolean value is False) 
    :PVInit: Initial value for PV (optional: default value is 0): used if PID_RT is ran first in the squence and no value of PV is available yet. 

    :method: discretisation method (optional: default value is 'EBD')
        EBD-EBD: EBD for integral action and EBD for derivative action
        EBD-TRAP: EBD for integral action and TRAP for derivative action
        TRAP-EBD: TRAP for integral action and EBD for derivative action
        TRAP-TRAP: TRAP for integral action and TRAP for derivative action

    The function "PID_RT" appends new values to the vectors "MV", "MVP", "MVI", and "MVD".
    The appended values are based on the PID algorithm, the controller mode, and feedforward.
    Note that saturation of "MV" within the limits [MVMin MVMax] is implemented with anti wind-up. 
    """

    if len(PV) == 0:
        E.append(SP[-1] - PVInit)
    else:
        E.append(SP[-1] - PV[-1]) 

    # Proportional action
    MVP.append(Kc*E[-1]) # à vérifier

    # Integral action
    if len(MVI) == 0:
        MVI.append((Kc*Ts/Ti)*E[-1]) # initialisation always with EBD
    else:
        if method == 'TRAP':
            MVI.append(MVI[-1] + (0.5*Kc*Ts/Ti)*(E[-1]+E[-2]))
        else:
            MVI.append(MVI[-1] + (Kc*Ts/Ti)*E[-1])

    # Derivative action # à vérifier
    Tfd = alpha*Td
    if len(MVD) == 0:
        #MVD.append(((Kc*Td)/(Tfd+Ts))*(E[-1]-E[-2])) # problem E[-2] n'existe pas
        MVD.append(((Kc*Td)/(Tfd+Ts))*E[-1]) # solution est de simplement supprimer E[-2]
    else:
        if method == 'TRAP':
            MVD.append(((Tfd-Ts*0.5)/(Tfd+Ts*0.5))*MVD[-1] + ((Kc*Td)/(Tfd+Ts*0.5))*(E[-1]-E[-2]))
        else:
            MVD.append((Tfd/(Tfd+Ts))*MVD[-1] + ((Kc*Td)/(Tfd+Ts))*(E[-1]-E[-2]))

    # Manual mode + anti wind-up 
    if Man[-1] == True:
        if ManFF:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1]
        else:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFF[-1]

##   Anti wind-up debug: Was implemented for a PI instead of a PID+FF, MVI did not decrease after MVFF > 0

    if (MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1]) > MVMax:
        MVI[-1] = MVMax - MVP[-1] - MVD[-1] - MVFF[-1]
    if (MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1]) < MVMin:
        MVI[-1] = MVMin - MVP[-1] - MVD[-1] - MVFF[-1]

    MV.append(MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1])
    
    # Saturation # à vérifier (peut-etre pas nécessaire)
    #MV[-1] = max(MVMin, min(MVMax, MV[-1]))
    

    return (MV, MVP, MVI, MVD)

def IMC(Kp, T1, T2, theta, gamma, order="FOPDT"):

    """
    :Kp: process gain
    :T1: First time constant [s]
    :T2: Second time constant [s] (Useful only if order="SOPDT")
    :theta: delay [s]
    :gamma: Regulator agressivity ([0.2-0.9], the smaller the more aggressive)
    :order: order of the model (optional: default value is 'FOPDT')
        FOPDT: First Order Plus Dead Time
        SOPDT: Second Order Plus Dead Time

    :return: IMC parameters (Kc, Ti, Td)
    """
    Tc = gamma*T1
    if(order=="FOPDT"):
        Kc = ((T1+(theta/2))/(Tc+(theta/2)))/Kp
        Ti = T1+(theta/2)
        Td = (T1*theta)/(2*T1+theta)
    elif(order=="SOPDT"):
        Kc = ((T1+T2)/(Tc+theta))/Kp
        Ti = T1+T2
        Td = (T1*T2)/(T1+T2)
    else:
        Kc = ((T1+(theta/2))/(Tc+(theta/2)))/Kp
        Ti = T1+(theta/2)
        Td = (T1*theta)/(2*T1+theta)
    return Kc, Ti, Td


class Controller:
    def __init__(self, parameters):
        self.parameters = parameters
        self.parameters['Kp'] = parameters['Kp'] if 'Kp' in parameters else 1.0
        self.parameters['Ti'] = parameters['Ti'] if 'Ti' in parameters else 0.0
        self.parameters['Td'] = parameters['Td'] if 'Td' in parameters else 0.0
        self.parameters['alpha'] = parameters['alpha'] if 'alpha' in parameters else 0.2



def margins(P, C, omega):

    """
    :P: Process object (e.g. SOPDT, FOPDT, etc.)
    :C: Controller object (e.g. PID, etc.)
    :omega: Frequency vector (e.g. np.logspace(-2, 2, 1000))

    :return: None. Opens a matplotlib figure with the Bode plot of the open loop transfer function L(s) = P(s)*C(s).
    """

    s = 1j*omega

    def find_nearest(array, value): # https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
        array = np.asarray(array)
        index = (np.abs(array - value)).argmin()
        result = array[index]
        return result, index

    def P_func(P, s): # work with "van der Grinten (SOPDT)" cause it is the best match
        Ptheta = np.exp(-P.parameters['theta']*s)
        PGain = P.parameters['Kp']*np.ones_like(Ptheta)
        PLag1 = 1/(P.parameters['Tlag1']*s + 1)
        PLag2 = 1/(P.parameters['Tlag2']*s + 1)
        #PLead1 = P.parameters['Tlead1']*s + 1
        #PLead2 = P.parameters['Tlead2']*s + 1
    
        Ps = np.multiply(Ptheta,PGain)
        Ps = np.multiply(Ps,PLag1)
        Ps = np.multiply(Ps,PLag2)
        #Ps = np.multiply(Ps,PLead1)
        #Ps = np.multiply(Ps,PLead2)
        return Ps

    def C_func(C, s): # PID -> MV = Kc * ( 1 + 1/(Ti*s) + (Td*s)/(alpha*Td*s + 1) ) * E
        CGain = C.parameters['Kc']
        CIntegral = 1 / (C.parameters['Ti'] * s)
        CDerivative = (C.parameters['Td'] * s) / (C.parameters['alpha'] * C.parameters['Td'] * s + 1)
    
        C_results = np.add(1, CIntegral, CDerivative)
        C_results = np.multiply(CGain, C_results)
        return C_results

    P_results = P_func(P, s)
    C_results = C_func(C, s)
    L_results = np.multiply(P_results, C_results)
    
    omega_c, omega_c_index = find_nearest(np.absolute(L_results), 1)
    omega_c = omega[omega_c_index]
    print(f"omega_c {omega_c} found at index {omega_c_index}")

    omega_u, omega_u_index = find_nearest(np.angle(L_results, True), -180)
    omega_u = omega[omega_u_index]
    print(f"omega_u {omega_u} found at index {omega_u_index}")

    fig, (ax_gain, ax_phase) = plt.subplots(2,1)
    fig.set_figheight(12)
    fig.set_figwidth(22)
    coord_delta = 1.2
    # Gain part
    om_u_y = 20*np.log10(np.abs(L_results[omega_u_index]))
    ax_gain.semilogx(omega,20*np.log10(np.abs(L_results)),label=r'$L(s)$')
    ax_gain.plot([-10, omega_c, omega_c], [0, 0, -1000], '-bo', label='omega_c')
    ax_gain.plot([-10, omega_u, omega_u], [om_u_y, om_u_y, -1000], '-go', label='omega_u')
    ax_gain.text(omega_c*coord_delta, 0, f'({omega_c:.5f}, {0})')
    ax_gain.text(omega_u*coord_delta, om_u_y, f'({omega_u:.5f}, {om_u_y:.5f})')
    gain_min = np.min(20*np.log10(np.abs(L_results)/5))
    gain_max = np.max(20*np.log10(np.abs(L_results)*5))
    ax_gain.set_xlim([np.min(omega), np.max(omega)])
    ax_gain.set_ylim([gain_min, gain_max])
    ax_gain.set_ylabel('Amplitude' + '\n $|L(j\omega)|$ [dB]')
    ax_gain.set_title('Bode plot of L = P*C')
    ax_gain.legend(loc='best')

    # Phase part
    om_c_y = (180/np.pi)*np.angle(L_results[omega_c_index])
    ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(L_results)),label=r'$L(s)$')
    ax_phase.plot([-10, omega_c, omega_c], [om_c_y, om_c_y, 0], '-bo', label='omega_c')
    ax_phase.plot([-10, omega_u, omega_u], [-180,-180,0], '-go', label='omega_u')
    ax_phase.text(omega_c*coord_delta, om_c_y, f'({omega_c:.5f}, {om_c_y:.5f})')
    ax_phase.text(omega_u*coord_delta, -180, f'({omega_u:.5f}, {-180})')
    ax_phase.set_xlim([np.min(omega), np.max(omega)])
    ph_min = np.min((180/np.pi)*np.unwrap(np.angle(L_results))) - 10
    ph_max = np.max((180/np.pi)*np.unwrap(np.angle(L_results))) + 10
    ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
    ax_phase.set_xlabel(r'Frequency $\omega$ [rad/s]')
    ax_phase.set_ylabel('Phase' + '\n $\,$'  + r'$\angle L(j\omega)$ [°]')
    ax_phase.legend(loc='best')

    return