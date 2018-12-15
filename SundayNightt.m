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
Q = 0.5;                              % volumetric flux [mÂ³/s]       % eg Neckar in TÃ¼bingen ~ 30 mÂ³/s
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
%§§§§§§§§§§§§§§§§§§§§ SPATIAL DISCRETIZATION  §§§§§§§§§§§§§§§§§§§§§§§§§§§§§
dx = 50;                % element length fix [m]
n = L_tot/dx ;        %  number of cells
x = dx/2:dx:L_tot-dx/2 ;  % central point of every cell
%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

% river geometry parameters in every cell
S_river = repmat(S_river,L_reach/dx, 1); S_river = S_river(:)';
B_0 = repmat(B_0,L_reach/dx,1); B_0 = B_0(:)';
s_bank = repmat(s_bank,L_reach/dx,1); s_bank = s_bank(:)';
          
%% CALCULATION of river geometry parameters
H = zeros(1,n);            % creates a vector for water depth [m]
A_c = zeros(1,n);          % creates a vector for cross sectional area [mÂ²]
A_s = zeros(1,n);          % creates a vector for surface area (water-atmosphere interface) [mÂ²]
q = zeros(1,n);            % creates a vector for specific discharge [m/s]

for i=1:n
    % Manning's equation solved for water depth H:    [QUAL2K maual eq(17)]
    if i == 1
        H(i)= (((Q*n_man)^(3/5))*(B_0(i)+2.*H_0*sqrt(s_bank(i).^2+1))^(2/5))/((S_river(i)^(3/10)).*(B_0(i)+s_bank(i).*H_0));
    else 
        H(i)= (((Q*n_man)^(3/5)).*(B_0(i)+2.*H(i)*sqrt(s_bank(i).^2+1)).^(2/5))/((S_river(i)^(3/10)).*(B_0(i)+s_bank(i).*H(i)));
    end    
A_c(i) = (B_0(i)+s_bank(i)*H(i))*H(i);  % cross-sectional area between elements [mÂ²] 
b_s(i) = H(i)/s_bank(i);                % length of water surface over sloped bank [m]
B_s(i) = B_0(i) + 2*b_s(i);             % total width of river at water surface [m]
A_s(i) = B_s(i)*L_reach;                % surface area (water-atmosphere interface) [mÂ²]
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
D = (0.011.*(q.^2).*(B_s.^2))./(Hmean.*Vsh);  % i think there is a time missing in the denominator? like this D has units of mÂ² but it should be mÂ²/s !?
Dnum = q*dx/2 ;  % numerical dispersion (QUAL2K.p18) 
Dm = zeros(1,length(D));
if Dnum <= D
    Dm = D - Dnum;
end
Dbulk= (Dm.*A_c)./dx  ;   % bulk diffusion coefficient

%§§§§§§§§§§§§§§§§§§§§ TIME DISCRETIZATION §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
dt = 60;                 % time step [s]
%te = 180*24*60*60;        % end time [s]
te=365*24*60*60;
t = 0:dt:te;
%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

%% EXTERNAL SOURCE PARAMETERS to be adjusted
pos = 120 ;        % cell where the flow is coming in
Q_s = 0.0001;      % external source discharge
T_s = 0;           % temperature of external source [K]

%% DATA IMPORT FOR RADIATION BALANCE
t_m=load('meteorological_input/t_m.txt');         % measurement times [d]  > we need measurement times in [s]
p_m=load('meteorological_input/p_m.txt');         % air pressure [kPa]  > we need air pres in [Pa] as input
T_a_m=load('meteorological_input/T_a_m.txt');     % air temperature[ï¿½C] > we need Temp in [K] as input
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
%t1_fix=800;
%t2_fix=250000;


t1_fix=173000; % for 90  days
t2_fix=259000; % for 180 days
t3_fix=389000; % for 270 days
t4_fix=514000; % for 356 days


figure (1)
plot (x/1000,Tw(t1_fix,:)-273.15,'r')
hold on
plot (x/1000,Tw(t2_fix,:)-273.15,'g')
plot (x/1000,Tw(t3_fix,:)-273.15,'b')
plot (x/1000,Tw(t4_fix,:)-273.15,'c')

% plot (x/1000,Tw_eq(150,:)-273.15,'b')
legend(sprintf('%6.2f days',t1_fix*dt/24/3600),sprintf('%6.2f days',t2_fix*dt/24/3600),sprintf('%6.2f days',t3_fix*dt/24/3600),sprintf('%6.2f days',t4_fix*dt/24/3600));
xlabel('x [Km]');
ylabel('T [{\circ}C]');
title('Space variation of Temperature at fixed time')
% ylim([265 288])
hold off
save RiverModel.mat;
save WaterTemperature Tw;

%% plot with fixed distance (year)
% x1_fix = 200;   % distance you want to show, in meters
% x2_fix = 30000; % second distance you want to show, in meters
x1_fix= 5000;  %T_w(column 100)
x2_fix= 15000; %T_w(column 300)
x3_fix= 25000; %T_w(column 500)
x4_fix= 30000; %T_w(column 600)

figure (2)

 plot (t/24/3600,Tw(:,100)-273.15,'r'), hold on
 plot (t/24/3600,Tw(:,300)-273.15,'g')
 plot (t/24/3600,Tw(:,500)-273.15,'b')
 plot (t/24/3600,Tw(:,600)-273.15,'k')
 
 plot (t/24/3600,T_a_m-273.15,'c')
% plot (Tw(1:1440,100)-273.15,'r'), 
% plot (t/24/3600,Tw_eq(:,150)-273.15,'b') 
 
 ylim([-8 40])
 xlim([0 350])
 
legend(sprintf('%6.2f m',x1_fix),sprintf('%6.2f m',x2_fix),sprintf('%6.2f m',x3_fix),sprintf('%6.2f m',x4_fix),'T_{air}')

xlabel('time [days]'); 
ylabel('T [{\circ}C]');
title('Time variation of Temperature at fixed position');
% datetick('x','mmm-dd');
% ylim([265 288])
hold off

%% Down sampling data

% here Tw data are downsampled 10 times. 
t_var = t/24/3600; 
t_var_downsampled = downsample(t_var,10);

y1 = Tw(:,100)-273.15; 
y1_downsampled = downsample(y1,10);

y2 = Tw(:,300)-273.15;
y2_downsampled = downsample(y2,10);

y3 = Tw(:,500)-273.15;
y3_downsampled = downsample(y3,10);

y4 = Tw(:,600)-273.15;
y4_downsampled = downsample(y4,10);

y5 = T_a_m-273.15;
y5_downsampled = downsample(y5,10);


% ----------------------- Curve fitting code ------------------------------

[xData1, yData1] = prepareCurveData( t_var_downsampled, y1_downsampled );
[xData2, yData2] = prepareCurveData( t_var_downsampled, y2_downsampled );
[xData3, yData3] = prepareCurveData( t_var_downsampled, y3_downsampled );
[xData4, yData4] = prepareCurveData( t_var_downsampled, y4_downsampled );
[xData5, yData5] = prepareCurveData( t_var_downsampled, y5_downsampled );


% Set up fittype and options.
ft = fittype( 'gauss4' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% Other Defaults:
%opts.Display = 'Off';
%opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];

% Fit model to data.
[fitresult1, gof1] = fit( xData1, yData1, ft, opts );
[fitresult2, gof2] = fit( xData2, yData2, ft, opts );
[fitresult3, gof3] = fit( xData3, yData3, ft, opts );
[fitresult4, gof4] = fit( xData4, yData4, ft, opts );
[fitresult5, gof5] = fit( xData5, yData5, ft, opts );

figure(3)

 plot (fitresult1,'r'), hold on
 plot (fitresult2,'g')
 plot (fitresult3,'b')
 plot (fitresult4,'k')
 plot (fitresult5,'c')

ylim([-8 40])
xlim([0 353])
 
xlabel('time [days]'); 
ylabel('T [{\circ}C]');
title('Time variation of Temperature at fixed position');
legend(sprintf('%6.2f m',x1_fix),sprintf('%6.2f m',x2_fix),sprintf('%6.2f m',x3_fix),sprintf('%6.2f m',x4_fix),'T_{air}')

hold off

%%


%$$$$$$$$$$$$$$$$$$$$$$$$$ my try to smooth the curve $$$$$$$$$$$$$$$$$$$$
% days=365;
% hoursperday= 24;
% coeff24hma= ones(1,hoursperday)/hoursperday;
% 
% avg24hTempC=filter(coeff24hma,1, tempC);
% 
% movingAvg_line=plot(days, avg24hTempC);






%% //////////////////// 3.OXYGEN MODEL ///////////////////////////////
%%

%  plot radiations for 3 month in detail


% figure(3)
% plot(t/24/3600,Tw(:,x1_fix/dx)-273.15,'r'), hold on






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
