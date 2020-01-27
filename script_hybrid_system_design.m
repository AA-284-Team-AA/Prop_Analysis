clear all
close all

%% Code for Hybrid System Design
% Based on G. Zilliac's slides on "Hybrid Propulsion System Design" 


%% INPUTS

% Step 1: Pick chamber pressure
p_c = 1.034e6; % Chamber pressure [Pa]
p_o = 101.5e3; % Atmospheric pressure [Pa]

% Step 2: Results from CEA 
of_ratio = 1.2879; % Optimal O/F ratio 
ISP = 208.63; %[s]
gamma = 1.4135;
c_star = 1661.781; %[m/s]

N = 1; % Number of ports
rho_f = 1180; % Fuel density [kg/m^3]
t_b = 20; %Burn time [s]

% Measured fuel constants for regression rate (Rabinovitch 2018)
% For use with SI units 
a = 8.96e-5; 
n = 0.35; 

ep = 2.5; % Nozzle area ratio, Ae/At

zeta_d = 1.07; % Discharge correction factor (ratio of actual to ideal mass flow rate)
zeta_v = 0.928; % Velocity correction factor (sqrt of energy conversion efficiency)

eff_c_star = 0.85; % c star efficiency

F = 50; % Initial thrust [N]

%OD = 3in
R_f = .0762/2; % Grain outer radius [m]

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

% Step 8: Calculates grain inner radius
syms R_i
temp1 = (R_f) ^(2*n + 1) - (R_i) ^ (2*n + 1);
temp2  = a * (2*n + 1) * (mdot_ox/pi) ^ n;
eqn = t_b == temp1 / temp2;
sol = vpasolve(eqn, R_i);

for j = 1: size(sol,1)
    if isreal(sol(j)) && sol(j)>0
        R_i = double(vpa(sol(j))) % Grain inner radius [m]
        break
    end
end

R_i = .3*.0254;
t = (0:.1:t_b);
for i = 1:size(t,2)
    R_inst(i) = (a*(2*n+1)*(mdot_ox/pi)^n*t(i) + R_i^(2*n+1))^(1/(2*n+1));
    final_R(i) = R_f; 
end
plot(t,R_inst)
hold on
plot(t,final_R);


% Calculates regression rate
G_ox = mdot_ox/(pi*(R_i^2)); %[kg/(m^2*s)]
rdot = (a * (G_ox)^n); %[m/s]

% Step 9: Calculates grain length
L = mdot_f / (2 * pi * N * R_i * rho_f * rdot) %[m]

% L = (10:.5:20).*.0254;
% R_i_solved = (mdot_f./(a.*(mdot_ox./pi).^n) .* 1./(rho_f.*2.*pi.*L)).^(1/(1-2*n));
% figure();
% plot(L./.0254,R_f - R_i_solved);


% step 10: Calculates required fuel and oxidizer mass
temp1 = pi * N * rho_f * L;
temp2 = (a * (2*n + 1) * (mdot_ox/(pi*N))^n * t_b + R_i^(2*n + 1))^(2/(2*n + 1));
m_f = temp1 * (temp2 - R_i^2); %[kg]
m_ox = mdot_ox * t_b; %[kg]
                    
% Calculates instantaneous O/F ratio
t = (0:.1:t_b);
inst_of_ratio = zeros(1,size(t,2));

for i = 1:size(t,2)
    temp1 = 1/(2*rho_f*L*a);
    temp2 = (mdot_ox/(pi*N))^(1-n);
    temp3 = (a * (2*n + 1) * (mdot_ox/(pi*N))^n * t(i) + R_i^(2*n + 1))^((2*n-1)/(2*n+1));
    inst_of_ratio(i) = temp1*temp2*temp3; % Instantaneous O/R ratio 
end






