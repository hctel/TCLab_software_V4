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
        TRAP: Trapezoïdal method
    
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

