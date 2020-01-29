import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel

def maintrade(fuel, oxi, testlims, N, P, testvar):
    ## inputs
    # range of chamber pressures to test (psi) (default: 400; can be list of numbers or range i.e. range(100,400,100))
    # testvar = "OF" # pick from OF, ER
    # oxi = "N2O"
    # fuel = "PMMA" # fuel is controlled by function input
    # OF = o/f ratio
    # ER = nozzle expansion ratio (Ae/At) (Default)
    # TESTLIMS: vector array (2 elements)--upper and lower limit of sweep in either OF or ER
    # TESTVAR: "OF" or "ER"

    # SETS UP and RUNS TRADES based on inputs
    pltname = oxi+fuel+"_"+testvar+"_" # prefix for plot image files
    if testvar == "ER":
        OFratio = testvar[1] # default=1.5
        supAR   = np.linspace(testlims[0],testlims[1],num=N,endpoint=True) # vary nozzle expansion ratio --this will make the supAR array
        ISP, Cstar, PCPE, cpcv = ARtrade(fuel,oxi,P,N,OFratio,supAR,pltname) # runs expansion ratio trade
        xval = supAR
    elif testvar == "OF":
        OFratio = np.linspace(testlims[0],testlims[1], num=N, endpoint=True) # vary O/F ratio (by mass) --this will make the OFratio array properly
        ISP, Cstar, PCPE, cpcv = OFtrade(fuel,oxi,P,N,OFratio,pltname) # runs O/F ratio trade
        xval = OFratio
    return ISP, Cstar, PCPE, cpcv, xval

# def findER(fu,ox,pcpe): # finds the ideal nozzle expansion ratio
#     C = CEA_Obj(oxName=ox, fuelName=fu,
#     isp_units='sec',
#     cstar_units='m/s') # define CEA object to operate on for rocketCEA
#     PCPE_fe = C.get_eps_at_PcOvPe(Pc=P,PcOvPe=pcpe)
#     return PCPE_fe

# the XXtrade functions work as templates for running any trade you might want. Just add more get_"" from rocketcea to work with more variables along with corresponding input
def ARtrade(fu,ox,P,N,OFratio,supAR,pltname): # expansion ratio trade
# fu: name of fuel (string, as defined in rocketcea documentation or newly added fuels)
# ox: name of ox (string, as defined in rocketcea documentation of newly added fuels)
# P: chamber pressure (either a single number, list of numbers, or range())
# N: number of desired supersonic area ratios (nozzle expansion ratio) to sweep over
# OFratio: fixed O/F ratio for this trade
# supAR: values of supersonic area ratios (length of this list must match value of N)
    C = CEA_Obj(oxName=ox, fuelName=fu,
    isp_units='sec',
    cstar_units='m/s') # define CEA object to operate on for rocketCEA
    if isinstance(P,int)==True: # if P is only one value
        y = 1
    else:
        y = len(P)
    # preallocate vars
    ISP     = np.zeros([y,supAR.shape[0]]) # isp
    Cstar   = np.zeros([y,supAR.shape[0]]) # cstar eff
    PCPE    = np.zeros([y,supAR.shape[0]]) # pc/pe
    cpcv    = np.zeros([y,supAR.shape[0]]) # ratio of specific heats in thrust chamber
    for x in range(y):
        if y==1:
            Pc = P # integers can't be called :(
            legends = str(Pc)
        else:
            Pc = P[x] # chamber pressure
            legends = P
        for i in range(N):
            ISP[x,i]    = C.get_Isp(Pc=Pc, MR=OFratio, eps=supAR[i]) # ISP vacuum
            Cstar[x,i]  = C.get_Cstar(Pc=Pc, MR=OFratio) # Cstar efficiency
            PCPE[x,i]   = C.get_PcOvPe(Pc=Pc, MR=OFratio, eps=supAR[i]) # Pc/Pe
            cpcv[x,i]   = C.get_Chamber_Cp(Pc=Pc, MR=OFratio, eps=supAR[i]) # cp/cv
    
    # generate plots for ISP, Cstar, and Pchamb/Pexit. Replace the last input with the vectory array of pressures
    # plots(supAR,ISP,"Ae/At","ISP (s)",  pltname+"isp.png"   , legends     ) # isp plot
    # plots(supAR,Cstar,"Ae/At","Cstar",  pltname+"cstar.png" , legends     ) # Cstar plot
    # plots(supAR,PCPE,"Ae/At","Pc/Pe",   pltname+"pcpe.png"  , legends     ) # Pc/Pe plot
    return ISP, Cstar, PCPE, cpcv

def OFtrade(fu,ox,P,N,OFratio,pltname): # O/F ratio trade (OFratio needs to be a vector array)
# fu: name of fuel (string, as defined in rocketcea documentation or newly added fuels)
# ox: name of ox (string, as defined in rocketcea documentation of newly added fuels)
# P: chamber pressure (either a single number, list of numbers, or range())
# N: number of desired O/F ratios to sweep over
# OFratio: values of O/F ratios (length of this list must match value of N)
# supAR: fixed nozzle expansion ratio
    C = CEA_Obj(oxName=ox, fuelName=fu,
    isp_units='sec',
    cstar_units='m/s') # define CEA object to operate on for rocketCEA
    if isinstance(P,int)==True: # if P is only one value
        y = 1
    else:
        y = len(P)
    # preallocate vars
    ISP     = np.zeros([y,OFratio.shape[0]]) # isp
    Cstar   = np.zeros([y,OFratio.shape[0]]) # cstar eff
    PCPE    = np.zeros([y,OFratio.shape[0]]) # pc/pe
    cpcv    = np.zeros([y,OFratio.shape[0]]) # ratio of specific heats in thrust chamber
    fe_pcpe = np.zeros([y,OFratio.shape[0]]) # nozzle area ratio for fully expanded flow
    for x in range(y):
        if y==1:
            Pc = P # integers can't be called :(
            legends = str(Pc)
        else:
            Pc = P[x] # chamber pressure
            legends = P
        pr = Pc/14.7 # pc/pe for fully expanded flo
        for i in range(N):
            fe_pcpe[x,:]= C.get_eps_at_PcOvPe(Pc=Pc, MR=OFratio[i], PcOvPe=pr)
            ISP[x,i]    = C.get_Isp(Pc=Pc, MR=OFratio[i], eps=fe_pcpe[x,i]) # ISP vacuum
            Cstar[x,i]  = C.get_Cstar(Pc=Pc, MR=OFratio[i]) # Cstar efficiency
            fe_pcpe[x,i]   = C.get_PcOvPe(Pc=Pc, MR=OFratio[i], eps=fe_pcpe[x,i]) # Pc/Pe
            cpcv[x,i]   = C.get_Chamber_Cp(Pc=Pc, MR=OFratio[i], eps=fe_pcpe[x,i]) # cp/cv
    
    # generate plots for ISP, Cstar, and Pchamb/Pexit
    # plots(OFratio,ISP,"O/F ratio","ISP (s)",    pltname+"isp.png"   , legends     ) # isp plot     
    # plots(OFratio,Cstar,"O/F ratio","Cstar",    pltname+"cstar.png" , legends     ) # Cstar plot        
    # plots(OFratio,PCPE,"O/F ratio","Pc/Pe",     pltname+"pcpe.png"  , legends     ) # Pc/Pe plot    
    return ISP, Cstar, fe_pcpe, cpcv

def plots(xvals,yvals,xname,yname,pltname,labels,plttit): # function to generate plots of the inputted variables
    plt.figure()
    if yvals.ndim==1:
        plt.plot(xvals,yvals[:], ms=10, label=str(labels))
        plt.xlim(min(xvals),max(xvals))
    else:
        for i in range(yvals.shape[0]): # can handle multiple lines (i.e. ISP vs. O/F at various chamber pressures)
            plt.plot(xvals,yvals[i,:], ms=10, label=str(labels[i]))
    plt.xlabel(xname)
    plt.ylabel(yname)
    plt.title(plttit)
    plt.legend(loc="lower right")
    plt.savefig(pltname)
    plt.close

### ANALYSIS SET UP ###
# defining fuel will add it's information to the master list of fuels to run rocketCEA with. 
# define PMMA
card_str = '''
fuel PMMA C 5 H 8 O 2
h,kj=-430.5 t(k)=299.82
''' # Greg Zilliac's recommendation for modeling PMMA
add_new_fuel('PMMA', card_str) # rocketCEA function to add PMMA to possible inputs

# define HTPB
card_str2 = '''
fuel HTPB C 7.3165 H 10.3360 O 0.1063
h,kj/mol= 456 t(k)=298.15 rho=0.9220
'''
add_new_fuel('HTPB', card_str2) # rocketCEA function to add HTPB to possible inputs

# define ABS (monomer of ABS)
card_str3 = '''
fuel ABS C 3 H 3 N 1
h,kj/mol=172 t(k)=299.82
'''
add_new_fuel('ABS (Monomer)', card_str3) # rocketCEA function to add HTPB to possible inputs

# define Paraffin
card_str4 = '''
fuel Paraffin C 32 H 66  wt%=100.00
h,kj/mol=-938 t(k)=298
'''
add_new_fuel('Paraffin', card_str4)


### BEGIN CODE TO RUN TRADES
testfs = ["PMMA","HTPB","ABS (Monomer)","Paraffin"]
testox = "GOX"
testvar = "OF" # pick OF or ER
N = 100 # number of points in trade study
P = 150 # chamber pressure, psi
pr = P/14.7 # pressure ratio Pc/Pe for fully expanded flow

ISPf    = np.zeros([len(testfs),N]) # ISP comparing all fuels
Cstarf  = np.zeros([len(testfs),N]) # Cstar comparing all fuels
PCPEf   = np.zeros([len(testfs),N]) # Pchamb/Pexit 
cpcvf   = np.zeros([len(testfs),N]) # ratio of specific heats
fe_pcpe = np.zeros([len(testfs),N]) # nozzle area ratio for fully expanded flow

pltlbls = [] # labels for plot legend
for i in range(len(testfs)): # labels for each line in plot
    pltlbls.append(testfs[i] + "/" + testox)
# currently setup for runs with only ONE chamber pressure selection
for f in range(len(testfs)):
    ISPf[f,:],Cstarf[f,:],PCPEf[f,:],cpcvf[f,:],xvar = maintrade(testfs[f],testox,[0.1,10],N,P,testvar) # currently overwriting xvar every time

# save plots of results
pltname = testox + "_Fuel Comparison_PC=" + str(P) + "_"
plttit = "P_c = " + str(P) + ",Fully Expanded"
if testvar == "OF":
    plots(xvar,ISPf,"O/F Ratio","ISP (s)",      pltname + "isp",pltlbls, "O/F Ratio vs. ISP, " + plttit)
    plots(xvar,Cstarf,"O/F Ratio","Cstar (m/s)",pltname + "cstar",pltlbls, "O/F Ratio vs. Cstar, " + plttit)
    plots(xvar,cpcvf,"O/F Ratio","Cp/Cv",       pltname + "cpcv",pltlbls, "O/F Ratio vs. Cp/Cv, " + plttit)
    plots(xvar,PCPEf,"O/F Ratio","Ae/At",       pltname + "aeat",pltlbls, "O/F Ratio vs. Ae/At, " + plttit)
elif testvar == "ER":
    plots(xvar,PCPEf,"Ae/At","Pc/Pe","Fuel Expansion Ratio Comparison.png",pltlbls)

# Note: for default case, Pc/Pe corresponding to fully expanded nozzle is 25.0
# why = np.reshape(np.array([[ISPf[0,:]],[Cstarf[0,:]]]),[2,50])
# plots(xvar,ISPf[0,:],"O/F","ISP and Cstar","ISP_Cstar Optimization",testfs)






 