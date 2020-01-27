clc
clear all
close all

%% Code for Hybrid System Design
% Based on G. Zilliac's slides on "Hybrid Propulsion System Design" 

%% Conversion factors

PA_TO_PSI = 0.000145038;
M_TO_IN = 39.3701;
KG_TO_LBS = 2.20462;

%% INPUTS

% Step 1: Pick chamber pressure
p_c = 150 / PA_TO_PSI; % Chamber pressure [Pa]
p_o = 101325; % Atmospheric pressure [Pa]

% Step 2: Results from CEA 
of_ratio = 1.2879; % Optimal O/F ratio 
ISP = 208.63; % [s]
gamma = 1.4135;
c_star = 1661.781; % [m/s]

N = 1; % Number of ports
rho_f = 1180; % Fuel density [kg/m^3]
t_b = 20; % Burn time [s]

% Measured fuel constants for regression rate (Rabinovitch 2018)
% For use with SI units 
a = 8.96e-5; 
n = 0.35; 

ep = 2.5; % Nozzle area ratio, Ae/At

zeta_d = 1.07; % Discharge correction factor (ratio of actual to ideal mass flow rate)
zeta_v = 0.928; % Velocity correction factor (sqrt of energy conversion efficiency)
eff_c_star = 0.85; % c star efficiency

F = 50; % Initial thrust [N]

R_f = 0.5 * 3/M_TO_IN; % Grain outer radius [m]

%% CALCULATIONS

% Step 3: Calculates nozzle exit pressure
syms p_e
temp1 = ((gamma + 1)/2)^(1/(gamma -1));
temp2 = (p_e/p_c)^(1/gamma);
temp3 = (gamma + 1)/(gamma - 1);
temp4 = 1 - (p_e/p_c)^((gamma -1)/gamma);
eqn = 0 == temp1*temp2*sqrt(temp3*temp4) - 1/ep;
[sol] = vpasolve(eqn,p_e);

for j = 1: size(sol,1)
    if isreal(sol(j)) && sol(j)>0
        p_e = double(vpa(sol(j))); % Nozzle exit pressure [Pa]
        break
    end
end

% Step 4: Calculates initial thrust coefficient 
temp1 = (2*gamma^2)/(gamma-1);
temp2 = (2/(gamma + 1))^((gamma+1)/(gamma-1));
temp3 = 1 - (p_e/p_c)^((gamma-1)/gamma);
temp4 = (p_e - p_o)*ep/p_c;
C_f = sqrt(temp1 * temp2 * temp3) - temp4; % Initial thrust coefficient

% Step 5: Efficiencies
eff_nozz = zeta_d * zeta_v; % Nozzle efficiency

% Step 6: Calculates initial throat area
A_t = F / (eff_nozz * C_f * p_c); % [m^2]

% Step 7: Calculates initial propellant mass flow rate 
mdot = p_c*A_t / (eff_c_star * c_star); % [kg/s]
mdot_f = mdot / (of_ratio + 1); % [kg/s]
mdot_ox = of_ratio * mdot_f; % [kg/s]

% R_f = .171/2*39.37;
% a = .223;
% n = .5;
% mdot_ox = 1.511*2.204;

% Step 8: Calculates maximum grain inner radius
syms R_i
temp1 = (R_f) ^(2*n + 1) - (R_i) ^ (2*n + 1);
temp2  = a * (2*n + 1) * (mdot_ox/pi) ^ n;
eqn = t_b == temp1 / temp2;
sol = vpasolve(eqn, R_i);

for j = 1: size(sol,1)
    if isreal(sol(j)) && sol(j)>0
        R_i_max = double(vpa(sol(j))); % Grain inner radius [m]
        break
    end
end

% R_i_max is to big because of the slow burning rate: pick R_i and check
% that it is less than R_i_max
R_i = 0.3 / M_TO_IN;
assert(R_i<R_i_max);

t = (0:.1:t_b);
for i = 1:size(t,2)
    R_i_inst(i) = (a*(2*n+1)*(mdot_ox/pi)^n*t(i) + R_i^(2*n+1))^(1/(2*n+1));
end

% Calculates initial regression rate
G_ox = mdot_ox/(pi*(R_i^2)); %[kg/(m^2*s)]
rdot = (a * (G_ox)^n); %[m/s]

% Step 9: Calculates grain length required to achieve initial thrust
L = mdot_f / (2 * pi * N * R_i * rho_f * rdot); %[m]
% or alternatively, obtain R_i given a range of grain length
L_range = (10:.5:20) / M_TO_IN;
R_i_range = ( mdot_f/(a*(mdot_ox/pi)^n) .* ...
            1./(rho_f*2*pi*L_range)) .^ (1/(1-2*n));

% step 10: Calculates required fuel and oxidizer mass
temp1 = pi * N * rho_f * L;
temp2 = (a * (2*n + 1) * (mdot_ox/(pi*N))^n * t_b + R_i^(2*n + 1))^(2/(2*n + 1));
m_f = temp1 * (temp2 - R_i^2); % [kg]
m_ox = mdot_ox * t_b; % [kg]
                    
% Calculates instantaneous O/F ratio
t = (0:.1:t_b);
inst_of_ratio = zeros(1,size(t,2));

for i = 1:size(t,2)
    temp1 = 1/(2*rho_f*L*a);
    temp2 = (mdot_ox/(pi*N))^(1-n);
    temp3 = (a * (2*n + 1) * (mdot_ox/(pi*N))^n * t(i) + R_i^(2*n + 1))^((2*n-1)/(2*n+1));
    inst_of_ratio(i) = temp1*temp2*temp3; % Instantaneous O/R ratio 
end

%% Plots
close all

% plot of R_i(t)
figure(1)
hold on
plot(t, R_i_inst * M_TO_IN)
plot(t, ones(length(t),1) * R_f * M_TO_IN, 'color', 'black')
xlabel('$t$ [s]'); ylabel('$R_i$ [in]')

% plot of R_i(L)
figure(2)
hold on 
plot(L_range * M_TO_IN, R_i_range * M_TO_IN)
xlabel('$L$ [in]'); ylabel('$R_i$ [in]')
xlim([L_range(1)*M_TO_IN; L_range(end)*M_TO_IN])















