# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel
## inputs
P = range(150,500,20) # range of chamber pressures to test

## add PMMA definition
card_str = '''
fuel PMMA C 5 H 8 O 2
h,kj=-430.5 t(k)=299.82
''' # Greg Zilliac's recommendation for modeling PMMA

# HTPB Definition (not used right now)
card_str2 = '''
fuel HTPB C 7.3165 H 10.3360 O 0.1063    wt%=100.00
h,cal= 1200.0 t(k)=298.15 rho=0.9220
'''

add_new_fuel('PMMA', card_str) # rocketCEA function to add PMMA to inputs


def rocket_vars(fu, ox, pcham):
    C = CEA_Obj(oxName=ox, fuelName=fu,
    isp_units='sec',
    cstar_units='m/s') # define CEA object to operate on for rocketCEA
    OFratio = np.linspace(0.1, 5., 50, endpoint=True) # OF ratio by mass
    ISP = np.zeros(OFratio.shape) # isp
    Cstar = np.zeros(OFratio.shape) # cstar efficiency
    PCPE = np.zeros(OFratio.shape) # p_chamber / p_exit 
    supAR   = 1   # supersonic area ratio (1 = converging nozzle only)
    for i in range(50):
        ISP[i] = C.get_Isp(pcham, MR=OFratio[i], eps=supAR) # ISP vacuum
        Cstar[i] = C.get_Cstar(pcham, MR=OFratio[i]) # Cstar efficiency
        PCPE[i] = C.get_Throat_PcOvPe(pcham, MR=OFratio[i]) # throat Pc/Pe (we will be using a converging nozzle)

    return ISP, Cstar, PCPE

pcham = 250
ISP_PMMA, Cstar_PMMA, PCPE_PMMA = rocket_vars('PMMA', 'GOX', pcham)
ISP_HTPB, Cstar_HTPB, PCPE_HTPB = rocket_vars('HTPB', 'GOX', pcham)

pltname = "GOXPMMA" + str(pcham) + ".png" #file name to save pressure plot
# parameters (unchanged)


# single parameter sensitivity study
N       = 50    # number of points
OFratio = np.linspace(0.1, 5., num=N, endpoint=True) # OF ratio by mass


fig, ax = plt.subplots(1,2)
ax[0].plot(OFratio, ISP_PMMA, color='green', ms=10)
ax[0].plot(OFratio, ISP_HTPB, ms=10)
ax[0].set_xlim(0., 4.)
ax[0].set_xlabel('OF ratio')
ax[0].set_ylabel('ISP [s]')

ax[1].plot(OFratio, Cstar_PMMA, color='green', ms=10)
ax[1].plot(OFratio, Cstar_HTPB, ms=10)
ax[1].set_xlim(0., 4.)
ax[1].set_xlabel('OF ratio')
ax[1].set_ylabel('$C^*$ [m/s]')

fig.set_size_inches(14., 7.)

plt.show()
