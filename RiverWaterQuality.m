% Case Studies 2018-19
% River Water Quality Model
% by Veronica, Alireza and Freya

%% To avoid conflicts: Who is working on which part?
% River geometry: Freya
% Radiation Balance: Alireza
% Oxygen model: Veronica


%%
%% MODEL
%%

%% ---------- 1. RIVER GEOMETRY ----------- 
%% COMMENTS on Geometry
% info: so far 3 reaches are computed, all have different hydraulic 
% parameters, this can be changed to more reaches if liked
% output parameters: (needed for further calculations) H, q, A_c, A_s with
% B_s (width @surface) for Dispersion

%% PARAMETERS to be adjusted
Q = 10;                              % volumetric flux [m³/s]       % eg Neckar in Tübingen ~ 30 m³/s
n = 0.03;                            % Manning roughness coefficient
S_river = [0.001 0.002 0.003];       % bottom slope [m/m]
%S_river = [0.001 0.002 0.004];         % bottom slope [m/m]
B_0 = [3 3 3];                       % bottom width [m]
%B_0 = [2.5 3 4];                       % bottom width [m] 
s_bank = [1 1 1];                    % slope of river banks [m/m] (dy/dy)
%s_bank = [1 2 0.5];                    % slope of river banks [m/m] (dy/dy)
H_0 = 2;                             % initial water depth in first reach
L_reach=10000;                       % reach length [m]

%% CALCULATION of water depth
H = zeros(3,1)';     % creates a vector for water depth [m]
A_c = zeros(3,1)';          % creates a vector for cross sectional area [m²]
A_s = zeros(3,1)';          % creates a vector for surface area (water-atmosphere interface) [m²]
q = zeros(3,1)';            % creates a vector for specific discharge [m/s]

for i=1:3
    % Manning's equation solved for water depth H:    [QUAL2K maual eq(17)]
    if i == 1
        H(i)= (((Q*n)^(3/5))*(B_0(i)+2.*H_0*sqrt(s_bank(i).^2+1))^(2/5))/((S_river(i)^(3/10)).*(B_0(i)+s_bank(i).*H_0));
    else 
        H(i)= (((Q*n)^(3/5)).*(B_0(i)+2.*H(i)*sqrt(s_bank(i).^2+1)).^(2/5))/((S_river(i)^(3/10)).*(B_0(i)+s_bank(i).*H(i)));
    end    
A_c(i) = (B_0(i)+s_bank(i)*H(i))*H(i);  % cross-sectional area between elements [m²] 
b_s(i) = H(i)/s_bank(i);                % length of water surface over sloped bank [m]
B_s(i) = B_0(i) + 2*b_s(i);             % total width of river at water surface [m]
A_s(i) = B_s(i)*L_reach;                % surface area (water-atmosphere interface) [m²]
q(i) = Q./A_c(i);                       % specific discharge [m/s]
end

%% -----------------------------------------------------
%% ------------ 2.TEMPERATURE MODEL --------------------
%% 2a. HEAT FLUXES WITHIN THE FLUIDS: ADVECTION, DISPERSION & CONDUCTION-------------------------

%% Parameters
% QUAL2K-Manual p.18
%info: in here the longitudinal dispersion (coefficient) is calculated for every reach
Hmean=H(1:3) ;                      % mean depth [m]  f
g=9.81 ;                            %[m/s^2]
Vsh=sqrt((g*Hmean.^2).*S_river);    % Shear velocity [m/s](QUAL2K.p18)
D=(0.011.*(q.^2).*(B_s.^2))./(q.*Vsh);

dx=(4.*D)./q   ;            %Courant number = 1, Neuman number =1/4
% dx =100 ;                 % element length fix [m]
dt = 4.*D./q.^2   ;         % with Courant number = 1 should be =dx/q;

Dbulk= (D.*A_c)./dx   ;

x=dx(1):dx(1):L_reach;
te = 60*60*2  ;
T_in= 300;
% T=T_in.*ones(1,length(x));
T = zeros(1,length(x));

  
 for t=0:dt(1):te
%========================= ADVECTION ====================================== 
T(2:end)=T(1:end-1);
T(1)= T_in;

%        figure
%       plot(x,T);
% %     hold on
%       ylim([0 320]);
% %     xlim([0 5000]);
% %     xlabel('x [m]');
% %     ylabel('T [K]');
%      title(sprintf('Temperature, t=%6.1f h',t/3600));
%       drawnow
    
%========================= DISPERSION ===================================== 
Hd=Dbulk(1).*(T(1:end-1)-T(2:end))./(A_c(1).*dx(1));
Hd =[0 Hd Hd(end)];                     % boundary conditions 
T = T+ dt(1).*(Hd(1:end-1)-Hd(2:end)); 

figure
plot(x,T);
ylim([0 320])
% %     hold on
%     xlabel('x [m]');
%     ylabel('T [K]');
title(sprintf('Temperature, t=%6.1f h',t/3600));
drawnow  
  end
 
 
% ========================= CONDUCTION =====================================
    


%% 2b. Sources & Sink Terms of Heat Balance -----------------------------

%% Exemplary measured Data
t_m   = [15; 46; 74; 105; 135; 166; 196; 227; 258; 288; 319; 349]*86400;                          % measurement times [s]
p_m   = [94540; 94380; 94980; 94430; 95280; 95150; 95330; 95110; 95490; 95470; 95670; 95480];     % air pressure [Pa]
T_a_m = [5.34; 1.23; 2.63; 7; 14.1; 16; 17.4; 16.9; 13.6; 10.7; 9; 3.2] + 273.15;                 % air temperature[K]
e_a_m = [530; 750; 880; 930; 1210; 1280; 1370; 1260; 1200; 940; 800; 767];                        % absolute humidty [Pa]
v_w_m = [1.04; 1.19; 2; 1.26; 1.67; 1; 2.88; 1.4; 1.5; 2.5; 1.3; 1.2];                            % wind speed[m/s]
B_m   = [0.79; 0.82; 0.65; 0.86; 0.76; 0.84; 0.51; 0.87; 0.98; 0.9; 0.75; 0.9];                   % cloud coverage [-]
H_G_m = [36.46; 63.66; 113.54; 125.5; 174.07; 201.22; 201.87; 183.68; 146.53; 89; 57.06; 55];     % global raditaion[W/m2]

%% PARAMETERS to be adjusted
% Characteristics of the lake
depth = 2;           % mean depth of the river [m]
depth_min = 1e-3;    % minimal depth to compute equilibrium depth
alpha = 0.05;        % albedo [-] value of 0.05 corresponds to ... surface (source:)

% Initial conditions
T_w0 = 13;              % initial river water temperature [°C]

tspan = [0:366]*86400; % Time span of one year in seconds


%% Model
% interpolation mode
mode='cubic';

% adjustment of input parameters
T_w0 = T_w0 + 273.15;   % initial river water temperature [K]

% Solve ODE
% ODE output is a time series of the water temperature
[t,T_w] = ode15s(@heatlakeode, tspan,T_w0,[],...
                  t_m,p_m,T_a_m,e_a_m,v_w_m,B_m,H_G_m,depth,alpha,mode);

% Convert [K] back to [°C]
T_w = T_w - 273.15;

%% Plot
figure
clf
plot(t/86400,T_w,'b');
title('Lake Temperature');
datetick('x','mmm-dd');
xlabel('Date');
ylabel('T_w [K]');
drawnow

% Compute with minimal depth => equilibrium temperature
% minimal depth was set to 1m < why?
[t,T_e] = ode15s(@heatlakeode, tspan,T_w0,[],...
               t_m,p_m,T_a_m,e_a_m,v_w_m,B_m,H_G_m,depth_min,alpha,mode);

% Convert [K] back to [°C]
T_e = T_e - 273.15;

hold on
plot(t/86400,T_e,':k',t/86400,...
    interp1(t_m,T_a_m,t,mode,'extrap')-273.15,'k--');
hold off
legend('T_{ini}=13K','T_{eq}','T_{air}','location','best') %,'T_{ini}=4K'
legend(gca,'boxoff')

%%
function dTdt = heatlakeode(t,T_w,t_m,p_m,T_a_m,e_a_m,v_w_m,B_m,H_G_m,...
                           depth,alpha,mode)
% ODE for the temperature model
% Constants
rho_w = 1000;       % mass density of water[kg/m^3]
W_ET = 2.5e6;       % specific enthalpy of volatilization of water [J/kg]
sigma = 5.67e-8;    % Stefan-Boltzmann-constant [W/(m^2 * K^-4)]
C_p = 4185;         % specific heat capacity of water [J/kg/K]

% interpolate the measured data to the current time point
p=interp1(t_m,p_m,t,'linear','extrap');
T_a=interp1(t_m,T_a_m,t,'linear','extrap');
e_a=interp1(t_m,e_a_m,t,mode,'extrap');
v_w=interp1(t_m,v_w_m,t,mode,'extrap');
B=interp1(t_m,B_m,t,mode,'extrap');
H_G=interp1(t_m,H_G_m,t,mode,'extrap');


%% Heat transfer in analogy to or coupled with mass tranfer
% Latent heat flux
theta = T_w - 273.15;                                    % water temperature in centigrades [�C]
e_sat = 6.111213*exp(17.5043*theta/(241.2 + theta))*100; % vapor pressure at saturation [Pa]
f = 5.44+2.19*v_w+0.24*(T_w-T_a);                        % wind factor [W/m2/Pa]
H_ET = f*(e_sat-e_a)/100;                                % latent heat flux [W/m2]

% Sensible heat flux
Bo = (T_w-T_a)/(e_sat-e_a)*p*6.5e-4; % Bowen ratio [-]
%H_c = H_ET*Bo;
kc = 0.203*sqrt(v_w);
H_c = kc*(T_w-T_a);

%% Heat transfer by radiation
% Heat budget due to short-wavelength radiation
H_insw = H_G ; % Measured radiation in the VIS spectrum (already accounting for cloud coverage)
H_outsw = alpha * H_insw; % reflected short-wave radiation

% Long wave radiation
H_lwout = sigma * T_w^4;                                 % heat radiation of lake
H_lwin = 6.8e-8*(e_a/100/T_a)^(1/7)*(1+0.17*B^2)*T_a^4;  % atmospheric backscattering

% Sum them up!
dTdt = (H_insw - H_outsw + H_lwin - H_lwout - H_ET - H_c)/(depth*rho_w*C_p); 
end




























