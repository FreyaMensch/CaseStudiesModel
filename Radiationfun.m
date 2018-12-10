function dTdt_R = Radiationfun(t,n,j,H,rho_w,C_p,T,alpha,sigma,p_m,T_a_m,R_H,v_w_m,H_G_m,B_m)
%RADIATIONFUN Summary of this function goes here
%   Detailed explanation goes here

for  i=1:n
theta(i) = T(i) - 273.15;                                    % water temperature in centigrades [°C]
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
 dTdt_R(i) = 1./(H(i)*rho_w*C_p).*(H_insw(i)- H_outsw(i)- H_lwout(i)- H_ET(i) - H_c(i)+ H_lwin(i));
end
end

