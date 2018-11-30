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

%% Load & pre-process measured meteorological data (exemplary)
% for ode input we need others
t_m=load('meteorological_input/t_m.txt');         % measurement times [d]  > we need measurement times in [s]
p_m=load('meteorological_input/p_m.txt');         % air pressure [kPa]  > we need air pres in [Pa] as input
T_a_m=load('meteorological_input/T_a_m.txt');     % air temperature[�C] > we need Temp in [K] as input
R_H=load('meteorological_input/RHumdity.txt');    % Relative humidty [%] > we need absolute humidity as input
v_w_m=load('meteorological_input/v_w_m.txt');     % wind speed[m/s]
H_G_m=load('meteorological_input/shortRadIn.txt');     % solar input radiation[w/m^2], Measured radiation in the VIS spectrum
B_m=(1-0.1).*rand(8724,1)+0.1;                    % Cloud coverage [-]

% meterological input data - daily mean
tm_n=floor(length(t_m));                % time in days [d] but as hourly values 
tm_n=(floor(length(t_m(4:end))/6));  % 6-hourly time as fractions of complete days [d]
t_m=(1:tm_n).*(60*60*24);            % time in seconds [s]
p_m = daily_mean(p_m) *1000;         % air pressure in [Pa]
T_a_m = daily_mean(T_a_m) + 273.15;  % air temperature in [K]
R_H=daily_mean(R_H);                 % Relative humidty [%]
v_w_m=daily_mean(v_w_m);
H_G_m=daily_mean(H_G_m);
B_m=daily_mean(B_m);                 % cloud coverage [-]

%% PARAMETERS to be adjusted
Q = 2;                              % volumetric flux [m³/s]       % eg Neckar in Tübingen ~ 30 m³/s
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
Vreachv(i)= A_c(i)*L_reach;                   % volume of reach [m³]
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

% Courant number = 1, Neuman number =1/4 :
% dx=(4.*D)./q   ;            
dx =10 ;                 % element length fix [m]
dt = min(4.*D./q.^2)   ;         % with Courant number = 1 should be =dx/q;
%dt = 60;                  % time fixed [s]

Dbulk= (D.*A_c)./dx   ;     % bulk diffusion coefficient

x=0:dx:L_reach ;
% make vectors for whole river from single reaches
A_c= repmat(A_c,[(length(x)-1) 1]);
A_c = A_c(:)';
H= repmat(H,[(length(x)-1) 1]);
H = H(:)';
q= repmat(q,[(length(x)-1) 1]);
q= q(:)';
D= repmat(D,[(length(x)-1) 1]);
D = D(:)';
Dnum=q*dx/4 ;  % numerical dispersion (QUAL2K.p18) /4 instead of /2
Dm=[0 0 0];
if Dnum<=D
    Dm=D-Dnum;
end

Dbulk= (Dm.*A_c)./dx  ;   % bulk diffusion coefficient

te = 2*60*60;               % end time [h]

t=0:dt:te;
T_in= 285;
n=L_tot/dx ;
% T=(T_in-10).*ones(1,length(x));
T = ones(1,n)*T_in;
T_eq=zeros(length(x),length(t));
 
Tw = ones(length(t),length(n));
dTdt_A = zeros(1,n);
%dTdt_R = zeros(length(V),1);

dt = 6*3600;    %for t=0:dt:te% dt has to be the nr of hours that the meteorological data is averaged over (6 hours in the present case)

for t=0:dt:te             
  for  i=2:n-1
%========================= ADVECTION ====================================== 
%  NOT RUNNING YET
%T(1)= T_in;
% T(2:end)=T(1:end-1);
dTdt_A(i)= (Q./A_c.*dx).*(T(i-1)-T(i));
  end
dTdt_A(1) = dTdt_A(1);     % boundary condition
T = T + dTdt*dt;
Tw(:,t) = T;
end

%========================= DISPERSION ===================================== 
% Hd=Dbulk(1).*(T(1:end-1)-T(2:end))./(A_c(1).*dx(1));
% Hd =[0 Hd Hd(end)];                     % boundary conditions 
% T = T+ dt(1).*(Hd(1:end-1)-Hd(2:end)); 
% 
% Tw(1:length(T))= T;



%========================= RADIATION ===================================== 
% I don't think subscripts need to be added as long as we dont want to save
% the variable in a matrix or as long as the variable is not changing over
% time

% %for  i=2:n-1                        % spacial iteration
% theta = T - 273.15;                                    % water temperature in centigrades [�C]
% e_sat = 6.111213.*exp(17.5043.*theta/(241.2 + theta)).*100; % vapor pressure at saturation [Pa] = Magnus equation
% e_a_m = (R_H(t)./100).*e_sat';
% 
% % Heat transfer in analogy to or coupled with mass tranfer
% % Latent heat flux
% f = 5.44+2.19*v_w_m(t)+0.24*(T-T_a_m(t));                        % wind factor [W/m2/Pa]
% H_ET = f*(e_sat-e_a_m(t))/100;                                % latent heat flux [W/m2]
% 
% % Sensible heat flux
% Bo = (T-T_a_m(t))/(e_sat-e_a_m(t))*p_m(t)*6.5e-4; % Bowen ratio [-]
% %H_c = H_ET*Bo;
% kc = 0.203*sqrt(v_w_m(t));
% H_c = kc*(T-T_a_m(t));
% 
% % Heat transfer by radiation
% % Heat budget due to short-wavelength radiation
% H_insw = H_G ; % Measured radiation in the VIS spectrum (already accounting for cloud coverage)
% H_outsw = alpha * H_insw; % reflected short-wave radiation
% 
% % Long wave radiation
% H_lwout = sigma * T^4;                                 % heat radiation of lake
% H_lwin = 6.8e-8*(e_a_m(t)/100/T_a_m(t))^(1/7)*(1+0.17*B_m(t)^2)*T_a_m(t)^4;  % atmospheric backscattering
% 
% % total temperature change over time = Sum of all heat fluxes
% dTdt = 1/(depth*rho_w*C_p)*(H_insw - H_outsw + H_lwin - H_lwout - H_ET - H_c);
% % dTdt = A_s/(depth* A_s*rho_w*C_p)*(H_insw - H_outsw + H_lwin - H_lwout - H_ET - H_c);


% %Plot
% figure(1)
% plot(x,T);
% ylim([270 290])
% % hold on
% xlabel('x [m]');
% ylabel('T [K]');
% title(sprintf('Temperature, t=%6.1f h',t/3600));
% drawnow

% collecting all values in a matrix:
% t=0:dt(1):te;
% T_eq(length(x),length(t))=T(x(i))

% end
 

%% 2b. Sources & Sink Terms of Heat Balance -----------------------------








%% PARAMETERS to be adjusted
% Characteristics of the lake
depth = H(1);           % mean depth of the river [m]
depth_min = 1e-3;    % minimal depth to compute equilibrium depth
alpha = 0.05;        % albedo [-] value of 0.05 corresponds to ... surface (source:)

% Initial conditions
 T_w0 = 1;                   % initial river water temperature [°C]
% %T_w0 = [4 6 8 12 16];
 T_w0 = T_w0 + 273.15;       % initial river water temperature [K]
%T_w0 = ;
tspan = [1:tm_n].* 86400;   % Time span of one year in seconds of 363days

%% Model
% interpolation mode
% mode='spline';                      % find reasonable mode.
% 
% % Solve ODE 
% % ODE output is a time series of the water temperature
% % for i= 1:5
% [t,T_w] = ode15s(@heatradiationODE, tspan,T_w0,[],...
%                   t_m,p_m,T_a_m,R_H,v_w_m,B_m,H_G_m,depth,alpha,mode);
% % end
%               
% % Compute with minimal depth => equilibrium temperature
% % minimal depth was set to 1m < why?
% [t,T_e] = ode15s(@heatradiationODE, tspan,T_w0,[],...
%                t_m,p_m,T_a_m,R_H,v_w_m,B_m,H_G_m,depth_min,alpha,mode);
% 
%           
% % Convert [K] back to [°C]
% T_e = T_e - 273.15;
% T_a_m = T_a_m - 273.15;  % air temperature in [C]
% T_w = T_w - 273.15;

%% Plot
figure
clf
plot(t/86400,T_w,'b');
datetick('x','mmm-dd');
hold on

% air temperature interpolation
plot(t/86400,interp1(t_m,T_a_m,t,mode,'extrap'),'k--');
%plot(t_m,T_a_m,'k--');
% plot(t2/86400,T_e,':k')
hold off

title('Water temperature in first reach (without advection & dispersion)');
xlabel('Date');
ylabel('T_w [C]');
legend('T_{w}=12C','T_{air}','location','best')    %,'T_{ini}=4K','T_{eq}'
legend(gca,'boxoff')

figure
plot(t/86400,H_G_m)
datetick('x','mmm-dd');
title('Incoming Radiation');
xlabel('Date');
ylabel('H [W/m²]');


%% RUBBISH / STUFF THAT MIGHT COME IN HANDY LATER AGAIN

% might be handy to produce a large matrix of certain extend
% % contour map of distribution
% cc=zeros(length(t),length(x));
% x=[0:0.01:2];
% for i=1:length(t)
%     for j=1:length(x)
%         cc(i,j)=1/sqrt(4*pi*D*t(i))*exp(-(x(j)-v*t(i))^2/(4*D*t(i)));
%     end
% end
% [X,Y]=meshgrid(x,t);


% might be come in handy later for faster reading of meteorological data:
% files = dir('textfiles/*.txt');
% for i=1:length(files)
%     eval(['load ' files(i).name ' -ascii']);
% end



