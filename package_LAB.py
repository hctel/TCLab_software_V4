import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

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

    if (MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1]) > MVMax:
        MVI[-1] = MVMax - MVP[-1] - MVD[-1] - MVFF[-1]
    if (MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1]) < MVMin:
        MVI[-1] = MVMin - MVP[-1] - MVD[-1] - MVFF[-1]

    MV.append(MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1])
    
    # Saturation # à vérifier (peut-etre pas nécessaire)
    #MV[-1] = max(MVMin, min(MVMax, MV[-1]))
    

    return (MV, MVP, MVI, MVD)

    '''
    # Input-output dynamics P(s)
    Delay_RT(MV,thetap,Ts,MVDelayp,MV0)
    FO_RT(MVDelayp,Kp,T1p,Ts,PV1p,0)
    FO_RT(PV1p,1,T2p,Ts,PV2p,0) 

    # Disturb  ce dynamics D(s)
    Delay_RT(DV - DV0*np.ones_like(DV),thetad,Ts,MVDelayd,0)
    FO_RT(MVDelayd,Kd,T1d,Ts,PV1d,0)
    FO_RT(PV1d,1,T2d,Ts,PV2d,0)

    PV.append(PV2p[-1] + PV2d[-1] + pV0-Kp*MV0)'
    '''

def IMC():
    #TODO: Implement
    pass