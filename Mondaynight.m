% Case Studies 2018-19
% River Water Quality Model
% by Veronica, Alireza and Freya

%% To avoid conflicts: Who is working on which part?
% River geometry: Freya
% Radiation Balance: Alireza
% Oxygen model: Veronica


%% Problems to solve:
% we have very different time scales.. in order to see dispersion we have to have very small time scales but on the year basis 

%%
%% MODEL
%%
clear all
close all

%% ---------- 1. RIVER GEOMETRY ----------- 
%% COMMENTS on Geometry
% info: so far 3 reaches are computed, all have different hydraulic 
% parameters, this can be changed to more reaches if liked
% output parameters: (needed for further calculations) H, q, A_c, A_s with
% B_s (width @surface) for Dispersion

%% Constants
rho_w = 1000;       % mass density of water[kg/m^3]
W_ET = 2.5e6;       % specific enthalpy of volatilization of water [J/kg]
sigma = 5.67e-8;    % Stefan-Boltzmann-constant [W/(m^2 * K^-4)]
C_p = 4185;         % specific heat capacity of water [J/kg/K]







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

%% PARAMETERS to be adjusted
Q = 0.1;                              % volumetric flux [m³/s]       % eg Neckar in Tübingen ~ 30 m³/s
n_man = 0.03;                            % Manning roughness coefficient
S_river = [0.001 0.002 0.003];       % bottom slope [m/m]
%S_river = [0.001 0.002 0.004];         % bottom slope [m/m]
B_0 = [4 5 6];                       % bottom width [m]
%B_0 = [2.5 3 4];                       % bottom width [m] 
s_bank = [1 1 1];                    % slope of river banks [m/m] (dy/dy)
%s_bank = [1 2 0.5];                    % slope of river banks [m/m] (dy/dy)
H_0 = 2;                             % initial water depth in first reach
L_tot = 30000;                       % total length of river [m]
reach_nr= 3;                         % number of reaches [-]
L_reach=L_tot/reach_nr;              % reach length [m]

% idea just for fun: S_river= 0.001*rand(10,1); B_m=(1-0.1).*rand(8724,1)+0.1;
% Tried, and it works > ask freya for code
%% CALCULATION of river geometry parameters
H = zeros(reach_nr,1)';            % creates a vector for water depth [m]
A_c = zeros(reach_nr,1)';          % creates a vector for cross sectional area [m²]
A_s = zeros(reach_nr,1)';          % creates a vector for surface area (water-atmosphere interface) [m²]
q = zeros(reach_nr,1)';            % creates a vector for specific discharge [m/s]

for i=1:reach_nr
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
Vreach(i)= A_c(i)*L_reach;                   % volume of reach [m³]
end
%% time it takes the water from start to end of river
%t_res = sum(L_reach./q);             % residence time of water parcel in river [seconds]
t_res = sum(L_reach./q)/3600;        % residence time of water parcel in river [hours]

%% -----------------------------------------------------
%% ------------ 2.TEMPERATURE MODEL --------------------
%% 2a. HEAT FLUXES WITHIN THE FLUIDS: ADVECTION, DISPERSION & CONDUCTION-------------------------

%% Parameters
% QUAL2K-Manual p.18
%info: in here the longitudinal dispersion (coefficient) is calculated for every reach
Hmean=H(1:reach_nr) ;                   % mean depth [m]
g=9.81 ;                                %[m/s^2]
Vsh=sqrt((g*Hmean.^2).*S_river);        % Shear velocity [m/s](QUAL2K.p18)
D=(0.011.*(q.^2).*(B_s.^2))./(Hmean.*Vsh);  % i think there is a time missing in the denominator? like this D has units of m² but it should be m²/s !?


dx=100;                % element length fix [m]
n = L_tot/dx ;     

x=dx/2:dx:L_reach-dx/2 ;
% make vectors for whole river from single reaches
A_c= repmat(A_c,[(length(x)) 1]);
A_c = A_c(:)';
H= repmat(H,[(length(x)) 1]);
H = H(:)';
q= repmat(q,[(length(x)) 1]);
q= q(:)';
D= repmat(D,[(length(x)) 1]);
D = D(:)';

Dnum=q*dx/2 ;  % numerical dispersion (QUAL2K.p18) 
Dm=zeros(1,length(D));
if Dnum<=D
    Dm=D-Dnum;
end

Dbulk= (Dm.*A_c)./dx  ;   % bulk diffusion coefficient

% Courant number = 1, Neuman number =1/4 :
% dx=(4.*D)./q   ;
% dt = min(4.*D./q.^2)   ;         % with Courant number = 1 should be =dx/q;

        
% dt=60;                % time fixed [s]
dt =120;    %for t=0:dt:te% dt has to be the nr of hours that the meteorological data is averaged over (6 hours in the present case)
% USING 6 HOURS THE ALGORITM FOR ADVECTION AND DISPERSION GOES CRAZY; WE
% HAVE TO US SMALLER TIME OR LARGER DX TO MAKE IT WORK.

te = 500*60*60;               % end time [s]
t=0:dt:te;

T_in = 284;
% n=L_tot/dx ;
T = ones(1,n)*272;
Tw = ones(length(t),n)*272;
 dTdt_A = zeros(1,n);
 dTdt_D = zeros(1,n);
 dTdt_S = zeros(1,n);
%dTdt_R = zeros(length(V),1);
Q_s= 0.001;    % external source discharge
T_s= 285;     % temperature of external source [K]
%% PARAMETERS to be adjusted
% Characteristics of the lake
% depth = H(1);           % mean depth of the river [m]
depth_min = 1e-3;    % minimal depth to compute equilibrium depth
alpha = 0.05;        % albedo [-] value of 0.05 corresponds to ... surface (source:)
%% Load & pre-process measured meteorological data (exemplary)
% for ode input we need others
t_m=load('meteorological_input/t_m.txt');         % measurement times [d]  > we need measurement times in [s]
p_m=load('meteorological_input/p_m.txt');         % air pressure [kPa]  > we need air pres in [Pa] as input
T_a_m=load('meteorological_input/T_a_m.txt');     % air temperature[�C] > we need Temp in [K] as input
R_H=load('meteorological_input/RHumdity.txt');    % Relative humidty [%] > we need absolute humidity as input
v_w_m=load('meteorological_input/v_w_m.txt');     % wind speed[m/s]
H_G_m=load('meteorological_input/shortRadIn.txt');     % solar input radiation[w/m^2], Measured radiation in the VIS spectrum
% B_m=(1-0.1).*rand(8724,1)+0.1;                    % Cloud coverage [-]
% interpolate to have the same time step as other processes (every 60 s)
[t_m, index] = unique(t_m); 
p_m= interp1(t_m*84600,p_m(index),t');
T_a_m= interp1(t_m*84600,T_a_m(index),t')+273.15;
R_H= interp1(t_m*84600,R_H(index),t');
v_w_m=interp1(t_m*84600,v_w_m(index),t');
H_G_m=interp1(t_m*84600,H_G_m(index),t');
B_m=(1-0.1).*rand(length(t),1)+0.1;                    % Cloud coverage [-]

%% 

theta = zeros(1,n);

for j=1:length(t)
    
%========================= ADVECTION ======================================
 for  i=2:n
dTdt_A(i)=(q(i)./dx).*(T(i-1)-T(i));
  end
dTdt_A(1)=(q(1)./dx).*(T_in-T(1));

%========================= DISPERSION =====================================
  for  i=2:n-1
dTdt_D(i)=(Dbulk(i-1)./(A_c(i)*dx)).*(T(i-1)-T(i)) + (Dbulk(i)./(A_c(i)*dx)).*(T(i+1)-T(i)) ;
  end
dTdt_D(1)=0 ;  %%%% ASK ABOUT THESE!!!
dTdt_D(n)=0 ;

%========================= EXTERNAL SOURCE ================================
for  i=100:100
W = rho_w*C_p*Q_s*T_s   ;

    dTdt_S(i)=W./(rho_w*C_p*A_c(i)*dx);
end
%========================= SUM OF ALL =====================================
for  i=1:n
theta(i) = T(i) - 273.15;                                    % water temperature in centigrades [�C]
e_sat(i) = 6.111213.*exp(17.5043.*theta(i)/(241.2 + theta(i))).*100; % vapor pressure at saturation [Pa] = Magnus equation
e_a_m(i) = (R_H(j)./100).*e_sat(i);

% Heat transfer in analogy to or coupled with mass tranfer
% Latent heat flux
f(i) = 5.44 + 2.19*v_w_m(j) + 0.24*(T(i) - T_a_m(j));                        % wind factor [W/m2/Pa]
H_ET(i) = f(i)*(e_sat(i)-e_a_m(i))/100;                      % latent heat flux [W/m2]

% Sensible heat flux
% H_c = H_ET*Bo;
kc = 0.203*sqrt(v_w_m(j));     % heat exchange coeff [W/(m^2 K)]
H_c(i) = kc.*(T(i) - T_a_m(j));

% Heat transfer by radiation
% Heat budget due to short-wavelength radiation
H_insw(i)= H_G_m(j) ; % Measured radiation in the VIS spectrum (already accounting for cloud coverage)
H_outsw(i) = alpha * H_insw(i); % reflected short-wave radiation

% % Long wave radiation
H_lwout(i) = sigma * T(i).^4;                                 % heat radiation
% vero(i)=6.8e-8.*(e_a_m(i)/100/T_a_m(j))^(1/7);
H_lwin(i) = 6.8e-8.*((e_a_m(i)/100/T_a_m(j))^(1/7))*(1+0.17*B_m(j)^2)*(T_a_m(j)^4);  % atmospheric backscattering
% 
% total temperature change over time = Sum of all heat fluxes
% dTdt(i) = 1./(H(i)*rho_w*C_p).*(H_insw(i) - H_outsw(i) + H_lwin(i) - H_lwout(i) - H_ET(i) - H_c(i));
 dTdt(i) = 1./(H(i)*rho_w*C_p).*(H_insw(i)- H_outsw(i)- H_lwout(i)- H_ET(i) - H_c(i)+ H_lwin(i));
end

T = T + dTdt_A*dt  + dTdt_D*dt + dTdt_S*dt + dTdt*dt;

Tw(j,:) = T;
  figure (1)
 plot (Tw(:,150))
ylim([265 288])





end


%========================= RADIATION ===================================== 
% I don't think subscripts need to be added as long as we dont want to save
% the variable in a matrix or as long as the variable is not changing over
% time

%% 


% 
% 
% % %Plot
% % figure(1)
% % plot(x,T);
% % ylim([270 290])
% % % hold on
% % xlabel('x [m]');
% % ylabel('T [K]');
% % hold on
% %  title(sprintf('Temperature, t=%6.1f h',t/3600));
% %  drawnow
% 
% % collecting all values in a matrix:
% % t=0:dt(1):te;
% % T_eq(length(x),length(t))=T(x(i))
% 
% end
% end
%  
% 
% %% Model
% % interpolation mode
% % mode='spline';                      % find reasonable mode.
% % 
%              
% % % Compute with minimal depth => equilibrium temperature
% % % minimal depth was set to 1m < why?
% % [t,T_e] = ode15s(@heatradiationODE, tspan,T_w0,[],...
% %                t_m,p_m,T_a_m,R_H,v_w_m,B_m,H_G_m,depth_min,alpha,mode);
% % 
% %           
% % % Convert [K] back to [°C]
% % T_e = T_e - 273.15;
% % T_a_m = T_a_m - 273.15;  % air temperature in [C]
% % T_w = T_w - 273.15;
% 
% %% Plot
% % figure
% % clf
% % plot(t/86400,T_w,'b');
% % datetick('x','mmm-dd');
% % hold on
% % 
% % % air temperature interpolation
% % plot(t/86400,interp1(t_m,T_a_m,t,mode,'extrap'),'k--');
% % %plot(t_m,T_a_m,'k--');
% % % plot(t2/86400,T_e,':k')
% % hold off
% % 
% % title('Water temperature in first reach (without advection & dispersion)');
% % xlabel('Date');
% % ylabel('T_w [C]');
% % legend('T_{w}=12C','T_{air}','location','best')    %,'T_{ini}=4K','T_{eq}'
% % legend(gca,'boxoff')
% % 
% % figure
% % plot(t/86400,H_G_m)
% % datetick('x','mmm-dd');
% % title('Incoming Radiation');
% % xlabel('Date');
% % ylabel('H [W/m²]');
% 
% 
% 
% 
% 
