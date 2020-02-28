import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
from scipy.optimize import fsolve


def get_exit_pressure(Pc, gamma, AR, guess=101325.):
    '''
    Find exit pressure from Pc, gam and AR.
    '''
    def f(x):
        temp1 = ((gamma + 1)/2) ** (1/(gamma -1))
        temp2 = (x/Pc) ** (1/gamma)
        temp3 = (gamma+1) / (gamma-1)
        temp4 = 1 - (x/Pc)**((gamma-1)/gamma);
        return temp1*temp2*np.sqrt(temp3*temp4) - 1/AR
    guess = 101325.
    Pe = fsolve(f, guess, xtol=1e-10)
    return Pe[0]


def get_CEA_results(C, Pc, Pe, OF, AR):
    PcOvPe = Pc/Pe
    _, gam_chamber = C.get_Chamber_MolWt_gamma(Pc=Pc, MR=OF, eps=AR)
    ISP = C.get_Isp(Pc=Pc, MR=OF, eps=AR)
    Cstar = C.get_Cstar(Pc=Pc, MR=OF)

    return gam_chamber, ISP, Cstar


def get_thrust_coeff(Pc, Pe, P0, gamma, AR):
    temp1 = (2*gamma**2) / (gamma-1)
    temp2 = (2/(gamma+1)) ** ((gamma+1)/(gamma-1))
    temp3 = 1 - (Pe/Pc) ** ((gamma-1)/gamma)
    temp4 = (Pe-P0)/Pc * AR
    return np.sqrt(temp1 * temp2 * temp3) + temp4


def get_mdots_from_Cstar(Pc, At, OF, Cstar, eff_Cstar=1.):
    mdot = Pc*At / (eff_Cstar*Cstar)
    mdot_fuel = mdot / (OF+1)
    mdot_ox = OF * mdot_fuel
    return mdot, mdot_ox, mdot_fuel


def get_radius_from_flow_rates(mdot_ox, mdot_fuel, rho_fuel, L, a, n):
    R = ( \
        mdot_fuel/(a*(mdot_ox/np.pi)**n) * 1./(rho_fuel*2*np.pi*L) \
        ) ** (1/(1-2*n))
    return R


def get_exit_velocity(Pc, Pe, at, gamma):
    Mach_exit = np.sqrt( \
        2/(gamma-1) * ( (Pc/Pe)**( (gamma-1)/gamma ) - 1) )
    tmp = (Pe/Pc)**( (gamma-1)/(2*gamma) ) * ((gamma+1)/2)**(1/2)
    return tmp * Mach_exit * at


def get_mdot_fuel(R, L, rho_fuel, mdot_ox, a, n):
    return 2*np.pi*R*L * rho_fuel * a*(mdot_ox/(np.pi*R**2))**n
