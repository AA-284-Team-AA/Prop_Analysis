# CEA_analysis

## Summary
Script to run fuel-ox analysis trades through RocketCEA. 

Goal is to facilitate quick runs and provide relevant plots.

## Instructions
Run FuelOxTrades.py for plot output. Below are the inputs to the script. All 4 possible fuels are shown below; number of tested fuels is unlimited (from 1 to infinity); number of chamber pressures per test is 1. Number of oxidizers per test is 1 maximum, 1 minimum.

```python
testfs = ["PMMA","HTPB","ABS (Monomer)","Paraffin"]
testox = "GOX"
testvar = "OF" # pick OF or ER
N = 100 # number of points in trade study
P = 150 # chamber pressure, psi
pr = P/14.7 # pressure ratio Pc/Pe for fully expanded flow
```

Code outputs plots with the name format "(Oxidizer)\_FuelComparison\_PC=(chamber pressure, in psi)\_(relevant variable).png". There are separate plots for ISP, Cstar, Cp/Cv, and expansion area ratio Ae/At. 

