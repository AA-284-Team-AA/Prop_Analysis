import numpy as np
import matplotlib.pyplot as plt
import pdb
import pickle
import utils
import configparser
import sys
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel


'''
Using Rabinovitch regression law, for L=18 [in]:
- nominal: OF[0] = 1.33 (to get R[0]=1 [in])
           Pc[0] = 215 [psi]
           throttling factor on mdot_ox: 4.5 to get to 20 [N]
Using Flora's regression law, for L=18 [in]
- nominal: OF[0] = 2.0 (to get R[0] of at least 0.15 [in])
           Pc[0] = 215 [psi]
           throttling factor on mdot_ox: 3.5 to get to 24 [N]
'''


# globals
FORMAT    = 'png'
PSI_TO_PA = 6894.7572931783
M_TO_IN   = 39.3701
G0        = 9.81
# define PMMA
card_str = '''
fuel PMMA C 5 H 8 O 2
h,kj=-430.5 t(k)=299.82
''' # Greg Zilliac's recommendation for modeling PMMA
add_new_fuel('PMMA', card_str) # rocketCEA function to add PMMA to possible inputs
CEA_object = CEA_Obj(oxName='GOX', fuelName='PMMA',
    isp_units='sec',
    cstar_units='m/s',
    sonic_velocity_units='m/s',
    temperature_units='k')


def get_initial_guess(config):
    Pc_0     = config.getfloat('Params', 'initial_pc') * PSI_TO_PA
    P0       = config.getfloat('Params', 'p0')
    OF_0     = config.getfloat('Params', 'initial_of')
    AR       = config.getfloat('Params', 'area_ratio')
    rho_fuel = config.getfloat('Params', 'rho_fuel')
    a        = config.getfloat('Params', 'a')
    n        = config.getfloat('Params', 'n')
    L        = config.getfloat('Params', 'fuel_length') / M_TO_IN
    Rt       = config.getfloat('Params', 'throat_radius')
    throat_area = np.pi * Rt**2

    gamma, ISP, Cstar = utils.get_CEA_results(CEA_object, Pc_0, P0, OF_0, AR)
    Pe = utils.get_exit_pressure(Pc_0, gamma, AR)
    CF = utils.get_thrust_coeff(Pc_0, Pe, P0, gamma, AR)
    mdot, mdot_ox, mdot_fuel = utils.get_mdots_from_Cstar(Pc_0,
        throat_area, OF_0, Cstar)
    R = utils.get_radius_from_flow_rates(mdot_ox, mdot_fuel, rho_fuel,
        L, a, n)
    return R, mdot_ox, Cstar


def initialize(config, R_0, mdot_ox_0, Cstar_0):
    n_timesteps = config.getint('Params', 'n_timesteps')
    burn_time   = config.getfloat('Params', 'burn_time')
    P0          = config.getfloat('Params', 'p0')
    rho_fuel    = config.getfloat('Params', 'rho_fuel')
    L           = config.getfloat('Params', 'fuel_length') / M_TO_IN
    a           = config.getfloat('Params', 'a')
    n           = config.getfloat('Params', 'n')
    Rt          = config.getfloat('Params', 'throat_radius')
    AR          = config.getfloat('Params', 'area_ratio')
    throat_area = np.pi * Rt**2

    t = np.linspace(0, burn_time, num=n_timesteps+1)

    R    = np.zeros(t.shape)
    R[0] = R_0

    mdot_ox    = np.zeros(t.shape)
    mdot_ox[0] = mdot_ox_0

    mdot_fuel    = np.zeros(t.shape)
    mdot_fuel[0] = utils.get_mdot_fuel(R[0], L, rho_fuel, mdot_ox[0],
        a, n)

    Pc    = np.zeros(t.shape)
    Pc[0] = Cstar_0 / throat_area * (mdot_ox[0]+mdot_fuel[0])

    OF    = np.zeros(t.shape)
    OF[0] = mdot_ox[0] / mdot_fuel[0]

    Cstar  = np.zeros(t.shape)
    thrust = np.zeros(t.shape)
    Pe     = np.zeros(t.shape)
    ISP    = np.zeros(t.shape)

    gamma, _, Cstar[0] = utils.get_CEA_results(CEA_object, Pc[0], P0,
        OF[0], AR)
    Pe[0] = utils.get_exit_pressure(Pc[0], gamma, AR)
    CF = utils.get_thrust_coeff(Pc[0], Pe[0], P0, gamma, AR)
    thrust[0] = Pc[0] * throat_area * CF
    ISP[0] = thrust[0] / G0 / (mdot_ox[0] + mdot_fuel[0])

    return t, R, mdot_ox, mdot_fuel, Pc, OF, Cstar, Pe, thrust, ISP


def get_desired_mdot_ox(config, tvec, mdot_ox):
    fac = config.getfloat('Params', 'throttling_factor')

    mdot_ox_0 = mdot_ox[0]
    mdot_ox_1, mdot_ox_2 = mdot_ox_0, mdot_ox_0/fac
    print(mdot_ox_1)
    print(mdot_ox_2)
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


def advance_R(a, n, R, mdot_ox, dt):
    Rdot = a * (mdot_ox / (np.pi * R**2)) ** n
    return R + dt*Rdot


def print_IC(t, R, mdot_ox, mdot_fuel, Pc, OF, Cstar, Pe, thrust):
    print('R[0]       = %f [in]' % (R[0]*M_TO_IN))
    print('mdot_ox[0] = %f [kg/s]' % mdot_ox[0])
    print('mdot_f[0]  = %f [kg/s]' % mdot_fuel[0])
    print('Pc[0]      = %f [psi]' % (Pc[0]/PSI_TO_PA))
    print('OF[0]      = %f' % OF[0])
    print('Cstar[0]   = %f' % Cstar[0])
    print('Pe[0]      = %f [Pa]' % Pe[0])
    print('thrust[0]  = %f [N]' % thrust[0])


def solve(config, t, R, mdot_ox, mdot_fuel, Pc, OF, Cstar, Pe, thrust, ISP):
    rho_fuel    = config.getfloat('Params', 'rho_fuel')
    L           = config.getfloat('Params', 'fuel_length') / M_TO_IN
    a           = config.getfloat('Params', 'a')
    n           = config.getfloat('Params', 'n')
    Rt          = config.getfloat('Params', 'throat_radius')
    AR          = config.getfloat('Params', 'area_ratio')
    P0          = config.getfloat('Params', 'p0')
    throat_area = np.pi * Rt**2

    N = t.shape[0]

    for i in range(N-1):
        dt = t[i+1] - t[i]
        R[i+1] = advance_R(a, n, R[i], mdot_ox[i], dt)

        mdot_fuel[i+1] = utils.get_mdot_fuel(R[i+1], L, rho_fuel,
            mdot_ox[i+1], a, n)

        mdot = mdot_fuel[i+1] + mdot_ox[i+1]
        Pc[i+1] = Cstar[i] / throat_area * mdot
        OF[i+1] = mdot_ox[i+1] / mdot_fuel[i+1]

        gamma, _, Cstar[i+1] = utils.get_CEA_results(CEA_object, Pc[i+1],
            P0, OF[i+1], AR)
        Pe[i+1] = utils.get_exit_pressure(Pc[i+1], gamma, AR, guess=Pe[i])
        CF = utils.get_thrust_coeff(Pc[i+1], Pe[i+1], P0, gamma, AR)
        thrust[i+1] = Pc[i+1]*throat_area*CF
        ISP[i+1] = thrust[i+1] / G0 / (mdot_ox[i+1] + mdot_fuel[i+1])


def save_results(t, R, mdot_ox, mdot_fuel, Pc, OF, Cstar, Pe, thrust, filename):
    res = {}
    res['t'] = t
    res['R'] = R
    res['mdot_ox'] = mdot_ox
    res['mdot_fuel'] = mdot_fuel
    res['Pc'] = Pc
    res['OF'] = OF
    res['Cstar'] = Cstar
    res['Pe'] = Pe
    res['thrust'] = thrust

    with open(filename + '.pkl', 'wb') as file:
        pickle.dump(res, file)


def plot_results(t, R, mdot_ox, mdot_fuel, Pc, thrust, ISP, filename):
    fig, ax = plt.subplots(3,3)
    ax = ax.flatten()
    for i in range(ax.shape[0]):
        ax[i].set_xlabel(r'$t$ [s]')
    i = 0

    ax[i].plot(t, mdot_ox)
    ax[i].set_ylabel(r'$\dot{m}_{ox}$ [${kg/m^3}$]')
    i += 1

    ax[i].plot(t, mdot_fuel)
    ax[i].set_ylabel(r'$\dot{m}_{fuel}$ [${kg/m^3}$]')
    i += 1

    ax[i].plot(t, mdot_ox/mdot_fuel)
    ax[i].set_ylabel(r'O/F')
    i += 1

    ax[i].plot(t, Pc / PSI_TO_PA)
    ax[i].set_ylabel(r'$P_c$ $[{psi}]$')
    i += 1

    ax[i].plot(t, thrust)
    ax[i].set_ylabel(r'$T$ [N]')
    i += 1

    ax[i].plot(t, 2*R*M_TO_IN)
    ax[i].set_ylabel(r'$D_{inner}$ [in]')
    i += 1

    ax[i].plot(t, Pe)
    ax[i].set_ylabel(r'$P_e$ [Pa]')
    i += 1

    ax[i].plot(t, mdot_ox / (np.pi * R**2))
    ax[i].set_ylabel(r'$G_{ox}$ [{kg/m^2/s}]')
    ax[i].set_ylim(0, 200.)
    i += 1

    ax[i].plot(t, ISP)
    ax[i].set_ylabel(r'$ISP$ [s]')
    i += 1

    fig.set_size_inches(24, 20)
    fig.savefig(filename+'.pdf', format='pdf')


def plot_results_separate(t, R, mdot_ox, mdot_fuel, Pc, thrust):
    fig, ax = plt.subplots()
    ax.plot(t, mdot_ox)
    ax.set_ylabel(r'$\dot{m}_{ox}$ [${kg/m^3}$]')
    fig.savefig('mdotox.'+FORMAT, format=FORMAT)

    fig, ax = plt.subplots()
    ax.plot(t, mdot_fuel)
    ax.set_ylabel(r'$\dot{m}_{fuel}$ [${kg/m^3}$]')
    fig.savefig('mdotfuel.'+FORMAT, format=FORMAT)

    fig, ax = plt.subplots()
    ax.plot(t, mdot_ox/mdot_fuel)
    ax.set_ylabel(r'O/F')
    fig.savefig('OF.'+FORMAT, format=FORMAT)

    fig, ax = plt.subplots()
    ax.plot(t, Pc / PSI_TO_PA)
    ax.set_ylabel(r'$P_c$ $[{psi}]$')
    fig.savefig('Pc.'+FORMAT, format=FORMAT)

    fig, ax = plt.subplots()
    ax.plot(t, thrust)
    ax.set_ylabel(r'$T$ [N]')
    fig.savefig('thrust.'+FORMAT, format=FORMAT)

    fig, ax = plt.subplots()
    ax.plot(t, 2*R*M_TO_IN)
    ax.set_ylabel(r'$D_{inner}$ [in]')
    fig.savefig('inner_diameter.'+FORMAT, format=FORMAT)

    fig, ax = plt.subplots()
    ax.plot(t, Pe)
    ax.set_ylabel(r'$P_e$ [Pa]')
    fig.savefig('Pe.'+FORMAT, format=FORMAT)

    fig, ax = plt.subplots()
    ax.plot(t, mdot_ox / (np.pi * R**2))
    ax.set_ylabel(r'$G_{ox}$ [{kg/m^2/s}]')
    fig.savefig('Gox.'+FORMAT, format=FORMAT)


def print_ranges(t, R, mdot_ox, mdot_fuel, Pc, thrust):
    print(30*'=')
    print('Printing ranges')
    print('Thrust: %f - %f' % (np.min(thrust), np.max(thrust)))
    print('mdot_ox: %f - %f' % (np.min(mdot_ox), np.max(mdot_ox)))
    print('Pc: %f - %f' % (np.min(Pc)/PSI_TO_PA, np.max(Pc)/PSI_TO_PA))
    print('O/F Ratio: %f - %f' % (np.min(mdot_ox/mdot_fuel),
        np.max(mdot_ox/mdot_fuel)))
    print(30*'=')


if __name__ == '__main__':
    filename = 'input.ini'
    if (len(sys.argv) == 2):
        filename = sys.argv[1]

    config = configparser.ConfigParser()
    config.read(filename)

    R_0, mdot_ox_0, Cstar_0 = get_initial_guess(config)
    t, R, mdot_ox, mdot_fuel, Pc, OF, Cstar, Pe, thrust, ISP = \
        initialize(config, R_0, mdot_ox_0, Cstar_0)
    get_desired_mdot_ox(config, t, mdot_ox)
    print_IC(t, R, mdot_ox, mdot_fuel, Pc, OF, Cstar, Pe, thrust)

    solve(config, t, R, mdot_ox, mdot_fuel, Pc, OF, Cstar, Pe, thrust, ISP)

    print_ranges(t, R, mdot_ox, mdot_fuel, Pc, thrust)


    if config.getboolean('Params', 'plot'):
        filename = config.get('Params', 'plot_file_name',
            fallback='time_analysis')
        # plot_results_separate(t, R, mdot_ox, mdot_fuel, Pc, thrust)
        plot_results(t, R, mdot_ox, mdot_fuel, Pc, thrust, ISP, filename)

    if config.getboolean('Params', 'save_results'):
        filename = config.get('Params', 'results_file_name',
            fallback='results')
        save_results(t, R, mdot_ox, mdot_fuel, Pc, OF, Cstar, Pe, thrust,
            filename)
