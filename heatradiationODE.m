function dTdt = heatradiationODE(t,T_w,t_m,p_m,T_a_m,R_H,v_w_m,B_m,H_G_m,...
                            depth,alpha,mode)
% ODE for the temperature model
% Constants
rho_w = 1000;       % mass density of water[kg/m^3]
W_ET = 2.5e6;       % specific enthalpy of volatilization of water [J/kg]
sigma = 5.67e-8;    % Stefan-Boltzmann-constant [W/(m^2 * K^-4)]
C_p = 4185;         % specific heat capacity of water [J/kg/K]

% calculate absolute humidity from relative humidity
theta = T_w - 273.15;                                    % water temperature in centigrades [�C]
e_sat = 6.111213.*exp(17.5043.*theta/(241.2 + theta)).*100; % vapor pressure at saturation [Pa] = Magnus equation
e_a_m = (R_H./100).*e_sat';

% interpolate the measured data to the current time point
p=interp1(t_m,p_m,t,'linear','extrap');
T_a=interp1(t_m,T_a_m,t,'linear','extrap');
e_a=interp1(t_m,e_a_m,t,mode,'extrap');
v_w=interp1(t_m,v_w_m,t,mode,'extrap');
B=interp1(t_m,B_m,t,mode,'extrap');
H_G=interp1(t_m,H_G_m,t,mode,'extrap');

%% Heat transfer in analogy to or coupled with mass tranfer
% Latent heat flux
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

% total temperature change over time = Sum of all heat fluxes
dTdt = 1/(depth*rho_w*C_p)*(H_insw - H_outsw + H_lwin - H_lwout - H_ET - H_c);
% dTdt = A_s/(depth* A_s*rho_w*C_p)*(H_insw - H_outsw + H_lwin - H_lwout - H_ET - H_c);
end

%% Plot
% figure
% plot(t_m,H_ET)
% hold on
% plot(t_m,H_c)
% plot(t_m,H_insw)
% plot(t_m,H_outsw)
% plot(t_m,H_lwin)
% plot(t_m,H_lwout)
% hold off
% title('Radiation Balance');
% xlabel('Date');
% ylabel('H [W/m²]');