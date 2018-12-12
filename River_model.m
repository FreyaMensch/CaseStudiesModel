% Case Studies 2018-19
% River Water Quality Model
% by Veronica, Alireza and Freya

%%
clear all
close all
%% /////////////////// 1. RIVER GEOMETRY //////////////////////////////// 
%% COMMENTS on Geometry
% info: so far 3 reaches are computed, all have different hydraulic 
% parameters, this can be changed to more reaches if liked
% output parameters: (needed for further calculations) H, q, A_c, A_s with
% B_s (width @surface) for Dispersion

%% PARAMETERS to be adjusted
Q = 0.5;                              % volumetric flux [m³/s]       % eg Neckar in Tübingen ~ 30 m³/s
n_man = 0.03;               % Manning roughness coefficient
L_tot = 30000;                       % total length of river [m]
reach_nr = 3;                         % number of reaches [-]
L_reach = L_tot/reach_nr; 
S_river = [0.001 0.002 0.003];       % bottom slope [m/m]
%S_river = [0.001 0.002 0.004];         % bottom slope [m/m]
B_0 = [5 6 7];                       % bottom width [m]
%B_0 = [2.5 3 4];                       % bottom width [m] 
s_bank = [1 1 1];                    % slope of river banks [m/m] (dy/dy)
%s_bank = [1 2 0.5];                    % slope of river banks [m/m] (dy/dy)
H_0 = 2;                             % initial water depth in first reach

% idea just for fun: S_river= 0.001*rand(10,1); B_m=(1-0.1).*rand(8724,1)+0.1;
% Tried, and it works > ask freya for code
%%
%�������������������� SPATIAL DISCRETIZATION  �����������������������������
dx = 50;                % element length fix [m]
n = L_tot/dx ;        %  number of cells
x = dx/2:dx:L_tot-dx/2 ;  % central point of every cell
%��������������������������������������������������������������������������

% river geometry parameters in every cell
S_river = repmat(S_river,L_reach/dx, 1); S_river = S_river(:)';
B_0 = repmat(B_0,L_reach/dx,1); B_0 = B_0(:)';
s_bank = repmat(s_bank,L_reach/dx,1); s_bank = s_bank(:)';
          
%% CALCULATION of river geometry parameters
H = zeros(1,n);            % creates a vector for water depth [m]
A_c = zeros(1,n);          % creates a vector for cross sectional area [m²]
A_s = zeros(1,n);          % creates a vector for surface area (water-atmosphere interface) [m²]
q = zeros(1,n);            % creates a vector for specific discharge [m/s]

for i=1:n
    % Manning's equation solved for water depth H:    [QUAL2K maual eq(17)]
    if i == 1
        H(i)= (((Q*n_man)^(3/5))*(B_0(i)+2.*H_0*sqrt(s_bank(i).^2+1))^(2/5))/((S_river(i)^(3/10)).*(B_0(i)+s_bank(i).*H_0));
    else 
        H(i)= (((Q*n_man)^(3/5)).*(B_0(i)+2.*H(i)*sqrt(s_bank(i).^2+1)).^(2/5))/((S_river(i)^(3/10)).*(B_0(i)+s_bank(i).*H(i)));
    end    
A_c(i) = (B_0(i)+s_bank(i)*H(i))*H(i);  % cross-sectional area between elements [m²] 
b_s(i) = H(i)/s_bank(i);                % length of water surface over sloped bank [m]
B_s(i) = B_0(i) + 2*b_s(i);             % total width of river at water surface [m]
A_s(i) = B_s(i)*L_reach;                % surface area (water-atmosphere interface) [m²]
q(i) = Q./A_c(i);                       % specific discharge [m/s]
end


%% time it takes the water from start to end of river
%t_res = sum(L_reach./q);             % residence time of water parcel in river [seconds]
t_res = sum(L_reach./q)/3600;        % residence time of water parcel in river [hours]

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
te = 60*24*60*60;        % end time [s]
t = 0:dt:te;
%��������������������������������������������������������������������������

%% EXTERNAL SOURCE PARAMETERS to be adjusted
pos = 120 ;        % cell where the flow is coming in
Q_s = 0.0001;      % external source discharge
T_s = 0;           % temperature of external source [K]

%% DATA IMPORT FOR RADIATION BALANCE
t_m=load('meteorological_input/t_m.txt');         % measurement times [d]  > we need measurement times in [s]
p_m=load('meteorological_input/p_m.txt');         % air pressure [kPa]  > we need air pres in [Pa] as input
T_a_m=load('meteorological_input/T_a_m.txt');     % air temperature[�C] > we need Temp in [K] as input
R_H=load('meteorological_input/RHumdity.txt');    % Relative humidty [%] > we need absolute humidity as input
v_w_m=load('meteorological_input/v_w_m.txt');     % wind speed[m/s]
H_G_m=load('meteorological_input/shortRadIn.txt');     % solar input radiation[w/m^2], Measured radiation in the VIS spectrum
% B_m=(1-0.1).*rand(8724,1)+0.1;                    % Cloud coverage [-]

% meterological input data - daily mean
% tm_n=floor(length(t_m));                % time in days [d] but as hourly values 
% tm_n=(floor(length(t_m(4:end))/6));  % 6-hourly time as fractions of complete days [d]
% t_m=(1:tm_n).*(60*60*24);            % time in seconds [s]
% p_m = daily_mean(p_m) *1000;         % air pressure in [Pa]
% T_a_m = daily_mean(T_a_m) + 273.15;  % air temperature in [K]
% R_H=daily_mean(R_H);                 % Relative humidty [%]
% v_w_m=daily_mean(v_w_m);
% H_G_m=daily_mean(H_G_m);
% B_m=daily_mean(B_m);                 % cloud coverage [-]


% interpolate to have the same time step as other processes (every dt)
[t_m, index] = unique(t_m); 
p_m= interp1(t_m*84600,p_m(index),t');
T_a_m= interp1(t_m*84600,T_a_m(index),t', 'spline')+273.15;
R_H= interp1(t_m*84600,R_H(index),t');
v_w_m=interp1(t_m*84600,v_w_m(index),t');
H_G_m=interp1(t_m*84600,H_G_m(index),t');
B_m=(1-0.1).*rand(length(t),1)+0.1;                    % Cloud coverage [-]


%% TIME LOOP TO CALCULATE TEMPERATURE (with Euler method)
T_in = 282;
T = ones(1,n)*282;   T_eq= ones(1,n)*282;
Tw = ones(length(t),n)*282; Tw_eq = ones(length(t),n)*282;
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
T_eq = T_eq + dTdt_A*dt + dTdt_D*dt + dTdt_S*dt ;  % T_eq is calculated by setting all air-water fluxes = 0
% 
Tw(j,:) = T;
Tw_eq(j,:) = T_eq;
end

%% plot with fixed time
t1_fix=800;      % index of Tw rows. Time in seconds will be t1_fix*dt
t2_fix=25000;    % index of Tw rows
figure (1)
plot (x/1000,Tw(t1_fix,:)-273.15,'r'), hold on
plot (x/1000,Tw(t2_fix,:)-273.15,'g')
% plot (x/1000,Tw_eq(150,:)-273.15,'b')
legend(sprintf('%6.2f days',t1_fix*dt/24/3600),sprintf('%6.2f days',t2_fix*dt/24/3600))
xlabel('x [Km]');
ylabel('T [{\circ}C]');
title('Space variation of Temperature at fixed time')
% ylim([265 288])
hold off

%% plot with fixed distance (year)
x1_fix = 1000;       % distance you want to show, in meters
x2_fix = 30000;       % second distance you want to show, in meters
figure (2)
plot (t/24/3600,Tw(:,x1_fix/dx)-273.15,'r'), hold on
plot (t/24/3600,Tw(:,x2_fix/dx)-273.15,'g')
plot(t/24/3600,T_a_m-273.15,'c')
legend(sprintf('%6.2f m',x1_fix),sprintf('%6.2f m',x2_fix),'T_{air}')
% plot (Tw(1:1440,100)-273.15,'r'), hold on

% plot (t/24/3600,Tw_eq(:,150)-273.15,'b')
xlabel('time [days]');
ylabel('T [{\circ}C]');
title('Time variation of Temperature at fixed position')
% datetick('x','mmm-dd');
% ylim([265 288])
hold off


%% //////////////////// 3.OXYGEN MODEL ///////////////////////////////
%%







%% EXTRA STUFF
% error = ones(1,n);       %[m] error matrix
% H     = ones(1,n);       %[m] initial water depth based on rough calculation
% while min(error) > 0.0001   %[m] minimum error for error matrix
% %     H_new =(((Q.*n_man).^(3/5)).*(B_0+2.*H.*sqrt(s_bank.^2+1)).^(2/5))./((S_river.^(3/10)).*(B_0+s_bank.*H));
%      H_new =(((Q.*n_man).^0.6).*(B_0+2.*H.*(s_bank.^2+1).^0.5).^0.4)./((S_river.^0.3).*(B_0+s_bank.*H));
%     error = abs((H_new-H)); %[m] absolute error between new H and old H
%     H = H_new;              %[m] update water depth matrix
% end
% 
% 
% Courant number = 1, Neuman number =1/4 :
% dx=(4.*D)./q   ;
% dt = min(4.*D./q.^2)   ;         % with Courant number = 1 should be =dx/q;

% meterological input data - daily mean
% tm_n=floor(length(t_m));                % time in days [d] but as hourly values 
% tm_n=(floor(length(t_m(4:end))/6));  % 6-hourly time as fractions of complete days [d]
% t_m=(1:tm_n).*(60*60*24);            % time in seconds [s]
% p_m = daily_mean(p_m) *1000;         % air pressure in [Pa]
% T_a_m = daily_mean(T_a_m) + 273.15;  % air temperature in [K]
% R_H=daily_mean(R_H);                 % Relative humidty [%]
% v_w_m=daily_mean(v_w_m);
% H_G_m=daily_mean(H_G_m);
% B_m=daily_mean(B_m);                 % cloud coverage [-]
