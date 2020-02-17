import numpy as np
import matplotlib.pyplot as plt
import pdb
import utils
from scipy.optimize import fsolve, excitingmixing, broyden1, newton_krylov, root
import matplotlib.pyplot as plt
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel


# globals
PSI_TO_PA = 6894.7572931783
M_TO_IN = 39.3701
# define PMMA
card_str = '''
fuel PMMA C 5 H 8 O 2
h,kj=-430.5 t(k)=299.82
''' # Greg Zilliac's recommendation for modeling PMMA
add_new_fuel('PMMA', card_str) # rocketCEA function to add PMMA to possible inputs
CEA_object = CEA_Obj(oxName='GOX', fuelName='PMMA',
    isp_units='sec',
    cstar_units='m/s',
    sonic_velocity_units='m/s')
REG_a = 8.96e-5
REG_n = 0.35
RHO_FUEL = 1180
ATM = 101325
EFF_NOZZLE = 1.0
EFF_CSTAR = 1.0

THROAT_AREA = 4.5861225022e-05
AREA_RATIO = 2.2949307022e+00
FUEL_LENGTH = 18 / M_TO_IN
BURN_TIME = 20
INITIAL_PC = 150*PSI_TO_PA
INITIAL_OF = 1.287
N_TIMESTEPS = 1000


def get_initial_mdot_ox():
    gamma, ISP, Cstar, at = utils.get_CEA_results(CEA_object, INITIAL_PC, ATM, \
        INITIAL_OF, AREA_RATIO)
    Pe = utils.get_exit_pressure(INITIAL_PC, gamma, AREA_RATIO)
    CF = utils.get_thrust_coeff(INITIAL_PC, Pe, ATM, gamma, AREA_RATIO)
    mdot, mdot_ox, mdot_fuel = utils.get_mdots_from_Cstar(INITIAL_PC,
        THROAT_AREA, INITIAL_OF, Cstar, eff_Cstar=EFF_CSTAR)
    R = utils.get_radius_from_flow_rates(mdot_ox, mdot_fuel, RHO_FUEL,
        FUEL_LENGTH, REG_a, REG_n)
    return R, mdot_ox, Cstar


def initialize(R_0, mdot_ox_0, Cstar_0):
    t = np.linspace(0, BURN_TIME, num=N_TIMESTEPS+1)

    R = np.zeros(t.shape)
    R[0] = R_0

    mdot_ox = np.zeros(t.shape)
    mdot_ox[0] = mdot_ox_0

    mdot_fuel = np.zeros(t.shape)
    mdot_fuel[0] = utils.get_mdot_fuel(R[0], FUEL_LENGTH, RHO_FUEL, mdot_ox[0],
        REG_a, REG_n)

    Pc = np.zeros(t.shape)
    Pc[0] = Cstar_0 / THROAT_AREA * (mdot_ox[0]+mdot_fuel[0])

    OF = np.zeros(t.shape)
    OF[0] = mdot_ox[0] / mdot_fuel[0]


    Cstar = np.zeros(t.shape)
    thrust = np.zeros(t.shape)

    gamma, _, Cstar[0], _ = utils.get_CEA_results(CEA_object, Pc[0], ATM,
        OF[0], AREA_RATIO)
    Pe = utils.get_exit_pressure(Pc[0], gamma, AREA_RATIO)
    CF = utils.get_thrust_coeff(Pc[0], Pe, ATM, gamma, AREA_RATIO)
    thrust[0] = Pc[0] * THROAT_AREA * CF

    return t, R, mdot_ox, mdot_fuel, Pc, OF, Cstar, thrust


def get_desired_mdot_ox(tvec, mdot_ox):
    mdot_ox_0 = mdot_ox[0]
    mdot_ox_1, mdot_ox_2 = mdot_ox_0, mdot_ox_0/2
    t1, t2 = 4, 15

    for i in range(t.shape[0]):
        if (t[i]<t1) or (t[i]>t2):
            mdot_ox[i] = mdot_ox_1
        elif (t[i]>=t1) and (t[i]<=t1+1):
            mdot_ox[i] = mdot_ox_1 + (mdot_ox_2-mdot_ox_1) * (t[i]-t1)
        elif (t[i]>=t2-1) and (t[i]<=t2):
            mdot_ox[i] = mdot_ox_2 + (mdot_ox_1-mdot_ox_2) * (t[i]-t2+1)
        else:
            mdot_ox[i] = mdot_ox_2


def advance_R(R, mdot_ox, dt):
    Rdot = REG_a * (mdot_ox / (np.pi * R**2)) ** REG_n
    return R + dt*Rdot


def solve(t, R, mdot_ox, mdot_fuel, Pc, OF, Cstar, thrust):
    N = t.shape[0]
    for i in range(N-1):
        dt = t[i+1] - t[i]
        R[i+1] = advance_R(R[i], mdot_ox[i], dt)

        mdot_fuel[i+1] = utils.get_mdot_fuel(R[i+1], FUEL_LENGTH, RHO_FUEL,
            mdot_ox[i+1], REG_a, REG_n)

        mdot = mdot_fuel[i+1] + mdot_ox[i+1]
        Pc[i+1] = Cstar[i] / THROAT_AREA * mdot
        OF[i+1] = mdot_ox[i+1] / mdot_fuel[i+1]

        gamma, _, Cstar[i+1], _ = utils.get_CEA_results(CEA_object, Pc[i+1],
            ATM, OF[i+1], AREA_RATIO)
        Pe = utils.get_exit_pressure(Pc[i+1], gamma, AREA_RATIO)
        CF = utils.get_thrust_coeff(Pc[i+1], Pe, ATM, gamma, AREA_RATIO)
        thrust[i+1] = Pc[i+1]*THROAT_AREA*CF


def plot_results(t, R, mdot_fuel, mdot_ox, Pc, thrust):
    fig, ax = plt.subplots(2,3)
    ax = ax.flatten()
    for i in range(ax.shape[0]):
        ax[i].set_xlabel(r'$t$ [s]')
    i = 0

    ax[i].plot(t, mdot_ox)
    ax[i].set_ylabel(r'$\dot{m}_\mathrm{ox}$ [$\SI{}{kg/m^3}$]')
    i += 1

    ax[i].plot(t, mdot_fuel)
    ax[i].set_ylabel(r'$\dot{m}_\mathrm{fuel}$ [$\SI{}{kg/m^3}$]')
    i += 1

    ax[i].plot(t, mdot_ox/mdot_fuel)
    ax[i].set_ylabel(r'O/F')
    i += 1

    ax[i].plot(t, Pc / PSI_TO_PA)
    ax[i].set_ylabel(r'$P_c$ $[\SI{}{psi}]$')
    i += 1

    ax[i].plot(t, thrust)
    ax[i].set_ylabel(r'$T$ [N]')
    i += 1

    ax[i].plot(t, 2*R*M_TO_IN)
    ax[i].set_ylabel(r'$D_\mathrm{inner}$ [in]')
    i += 1


if __name__ == '__main__':
    R_0, mdot_ox_0, Cstar_0 = get_initial_mdot_ox()
    t, R, mdot_ox, mdot_fuel, Pc, OF, Cstar, thrust = \
        initialize(R_0, mdot_ox_0, Cstar_0)
    get_desired_mdot_ox(t, mdot_ox)
    solve(t, R, mdot_ox, mdot_fuel, Pc, OF, Cstar, thrust)
    plot_results(t, R, mdot_fuel, mdot_ox, Pc, thrust)
    plt.show()
