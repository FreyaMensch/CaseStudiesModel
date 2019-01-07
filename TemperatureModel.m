%% //////////////////// 2.TEMPERATURE MODEL ///////////////////////////////
%% Constant Parameters
rho_w = 1000;       % mass density of water[kg/m^3]
W_ET = 2.5e6;       % specific enthalpy of volatilization of water [J/kg]
sigma = 5.67e-8;    % Stefan-Boltzmann-constant [W/(m^2 * K^-4)]
C_p = 4185;         % specific heat capacity of water [J/kg/K]
depth_min = 1e-3*ones(n,1);    % minimal depth to compute equilibrium depth
alpha = 0.05;        % albedo [-] value of 0.05 corresponds to ... surface (source:)

% Longitudinal dispersion (coefficient) is calculated for every cell
Hmean = H ;                   % mean depth [m]
g = 9.81 ;                                % gravity [m/s^2]
Vsh = sqrt((g*Hmean.^2).*S_river);        % Shear velocity [m/s](QUAL2K.p18)
D = (0.011.*(q.^2).*(B_s.^2))./(Hmean.*Vsh);  % i think there is a time missing in the denominator? like this D has units of m² but it should be m²/s !?
Dnum = q*dx/2 ;  % numerical dispersion (QUAL2K.p18) 
Dm = zeros(1,length(D));
if Dnum <= D
    Dm = D - Dnum;
end
Dbulk= (Dm.*A_c)./dx  ;   % bulk diffusion coefficient

%�������������������� TIME DISCRETIZATION ���������������������������������
dt = 60;                 % time step [s]
te = 360*24*60*60;        % end time [s]
t = 0:dt:te;
%��������������������������������������������������������������������������

%% EXTERNAL SOURCE PARAMETERS to be adjusted
pos = 120 ;        % cell where the flow is coming in
Q_s = 0.001;      % external source discharge
T_s = 285;           % temperature of external source [K]

%% DATA IMPORT FOR RADIATION BALANCE
t_m=load('meteorological_input/t_m.txt')*(24*3600);         % measurement times [d]  > we need measurement times in [s]
p_m=load('meteorological_input/p_m.txt')*1000;         % air pressure converted from [kPa] to [Pa]
T_a_m=load('meteorological_input/T_a_m.txt')+ 273.15;  % air temperature converted from [°C] to [K]
R_H=load('meteorological_input/RHumdity.txt');    % Relative humidty [%] > we need absolute humidity as input
v_w_m=load('meteorological_input/v_w_m.txt');     % wind speed[m/s]
H_G_m=load('meteorological_input/shortRadIn.txt');     % solar input radiation[w/m^2], Measured radiation in the VIS spectrum
B_m=(1-0.1).*rand(length(t),1)+0.1;                    % Cloud coverage [-]


% interpolate to have the same time step as other processes (every dt)
[t_m, index] = unique(t_m); 
p_m= interp1(t_m,p_m(index),t');
T_a_m= interp1(t_m,T_a_m(index),t', 'spline');
R_H= interp1(t_m,R_H(index),t');
v_w_m=interp1(t_m,v_w_m(index),t');
H_G_m=interp1(t_m,H_G_m(index),t');
B_m=interp1(t_m,B_m(index),t');


%% TIME LOOP TO CALCULATE TEMPERATURE (with Euler method)
T_in = 282;
T = ones(1,n)*282;   %T_eq= ones(1,n)*282;
Tw = ones(length(t),n)*282; %Tw_eq = ones(length(t),n)*282;
dTdt_A = zeros(1,n); dTdt_D = zeros(1,n); dTdt_S = zeros(1,n);

for j=1:length(t)
%========================= ADVECTION ======================================
dTdt_A = Advectionfun(n,q,dx,T,T_in);
%========================= DISPERSION =====================================
dTdt_D = Dispersionfun(n,Dbulk,A_c,dx,T);
%========================= EXTERNAL SOURCE ================================
for  i=pos:pos             % point source (inserted in only one cell)
     W = rho_w*C_p*Q_s*T_s   ;
    dTdt_S(i)=W./(rho_w*C_p*A_c(i)*dx);
end
%========================= RADIATION ======================================
dTdt_R = Radiationfun(t,n,j,H,rho_w,C_p,T,alpha,sigma,p_m,T_a_m,R_H,v_w_m,H_G_m,B_m);
% dTdt_R_eq = Radiationfun(t,n,j,depth_min,rho_w,C_p,T,alpha,sigma,p_m,T_a_m,R_H,v_w_m,H_G_m,B_m);
%========================= SUM OF ALL =====================================
T = T + dTdt_A*dt + dTdt_D*dt + dTdt_S*dt + dTdt_R*dt;
%T_eq = T_eq + dTdt_A*dt + dTdt_D*dt + dTdt_S*dt ;  % T_eq is calculated by setting all air-water fluxes = 0
% 
Tw(j,:) = T;
%Tw_eq(j,:) = T_eq;
end
