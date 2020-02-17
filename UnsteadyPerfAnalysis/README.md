# Temporal performance analysis

## Step 1: initialization

From desired initial chamber pressure and OF ratio, the code computes $R_i^0$, $\dot{m}_\mathrm{ox}^0$ and $C^*_\mathrm{guess}$. The code does not care about outer radius and the user should check that it is always greater than the inner radius.

## Step 2: solve

Then, the input for the temporal analysis is a profile $\dot{m}_\mathrm{ox}^i$. The sequence of operations are:

1. Get $R_i^{i+1}$ using a basic forward Euler
2. Get $\dot{R}_i^{i+1}$, $\dot{m}_\mathrm{fuel}^{i+1}$ and $\mathrm{O/F}^{i+1}$
3. Get $p_c^{i+1} = C^{*,i} \dot{m}^{i+1} / A_t$ (lagged $C^*$)
4. Recompute $\gamma$ and $C^{*,i+1}$ from $p_c^{i+1}$, $\mathrm{O/F}^{i+1}$
5. Get $C_F$ from $p_c^{i+1}$ and $\gamma$
6. Get $T^{i+1} = p_c^{i+1} A_t C_F$

$C^*_\mathrm{guess}$ is used for the inialization following the same procedure.