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
        H(i)= (((Q*n)^(3/5))*(B_0(i)+2.*H_0*sqrt(s_bank(i).^2+1))^(2/5))/((S_river(i)^(3/10)).*(B_0(i)+s_bank(i).*H_0));
    else 
        H(i)= (((Q*n)^(3/5)).*(B_0(i)+2.*H(i)*sqrt(s_bank(i).^2+1)).^(2/5))/((S_river(i)^(3/10)).*(B_0(i)+s_bank(i).*H(i)));
    end    
A_c(i) = (B_0(i)+s_bank(i)*H(i))*H(i);  % cross-sectional area between elements [m²] 
b_s(i) = H(i)/s_bank(i);                % length of water surface over sloped bank [m]
B_s(i) = B_0(i) + 2*b_s(i);             % total width of river at water surface [m]
A_s(i) = B_s(i)*L_reach;                % surface area (water-atmosphere interface) [m²]
q(i) = Q./A_c(i);                       % specific discharge [m/s]
V(i)= A_c(i)*L_reach;                   % volume of reach [m³]
end

    
%% time it takes the water from start to end of river
%t_res = sum(L_reach./q);             % residence time of water parcel in river [seconds]
t_res = sum(L_reach./q)/3600;        % residence time of water parcel in river [hours]
%% -----------------------------------------------------
%% ------------ 2.TEMPERATURE MODEL --------------------
%% 2a. HEAT FLUXES WITHIN THE FLUIDS: ADVECTION, DISPERSION & -----------------------

%% Parameters
% QUAL2K-Manual p.18
%info: in here the longitudinal dispersion (coefficient) is calculated for every reach
Hmean=H(1:reach_nr) ;                   % mean depth [m]
g=9.81 ;                                %[m/s^2]
Vsh=sqrt((g*Hmean.^2).*S_river);        % Shear velocity [m/s](QUAL2K.p18)
D=(0.011.*(q.^2).*(B_s.^2))./(Hmean.*Vsh);  % 

% Courant number = 1, Neuman number =1/2 :
dx=(4.*D)./q  ;            
dt = 4.*D./q.^2  ;   % with Courant number = 1 should be =dx/q;
dx = round(min(dx));
dt = round(min(dt));
te = 6.3*60*60;               % end time [h]
tspan=0:dt:te;
t=0:dt:te;

%   dx = 100 ;
%  dt = 60 ;

x=0:dx:L_reach ;
A_c=[A_c(1)*ones(1,length(x)),A_c(2)*ones(1,length(x)), A_c(3)*ones(1,length(x))];
H=[H(1)*ones(1,length(x)),H(2)*ones(1,length(x)), H(3)*ones(1,length(x))];
q=[q(1)*ones(1,length(x)),q(2)*ones(1,length(x)), q(3)*ones(1,length(x))];
D=[D(1)*ones(1,length(x)),D(2)*ones(1,length(x)), D(3)*ones(1,length(x))];
%

Dnum=q*dx/4 ;  % numerical dispersion (QUAL2K.p18) /4 instead of /2
Dm=[0 0 0];
if Dnum<=D
    Dm=D-Dnum;
end

Dbulk= (Dm.*A_c)./dx  ;   % bulk diffusion coefficient
V= A_c.*dx ;                 % volume of element 
T_0 = [283 273.*ones(1, length(V)-1)];

%========================= ODE SOLVE ====================================== 

    
[t,T]=ode45(@AdvDisFun,tspan,T_0,[],x,Q,V,Dbulk) ;

T=T-273.15 ;
figure

plot([0:dx:L_tot]/1000,T(end,:));
% hold on
% plot([0:dx:L_tot]/1000,T(2500,:));
% plot([0:dx:L_tot]/1000,T(4500,:));
% xlabel('x [Km]');
% ylabel('T [�C]');
% hold on

figure(2)
 plot(t/3600,T);
%  hold on
%  plot(t/3600,T(:,500));
%   plot(t/3600,T(:,1600));
 xlabel('t [h]');
ylabel('T [�C]');
 % figure(2)
% plot(t,T{k}(:,50));

%