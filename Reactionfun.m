function [dcdt_react_O2,dcdt_react_bod] = Reactionfun(~,j,n,H_G_m,c_x,c_x_bod,H,Dm,q,Tw,csat)
%REACTIONFUN Summary of this function goes here
%   Detailed explanation goes here
% rate coefficient for photosynthesis
k_photo = 0.05/86400; % g/(Ws) (problem 16)
% maximum respiration rate
k_resp_max = 200/86400; %g/m^2/s  (problem 16)
% Michaelis-Menten coefficients
KO2 = 0.5;
KBOD = 2;
dcdt_react_O2 = zeros(1,n);
dcdt_react_bod = zeros(1,n);

     theta=1.056; % temperature coefficient. This has a value of 1.056 in general
%     k_20=1.104; % BOD rate constant determined at 20oC, day-1
%     k_e=0.1;    % the light extinction coefficient, [1/m],, from the paper(light coeff coef)page.24 i just put it in the latex leterature,but not mentioned in the text yet!
%     K_lp= 0.2 ;    % bottom algae parameter [w/m2] !!!
%     r_max= ;   % the maximum photosynthesis rate,[g/m^2/s]
    
for  i=1:n
% Oxygen input via PHOTOSYNTHESIS
r_photo = *k_photo*H_G_m(j)./H;  % (problem 16) has to be made Temperature dependant!!

%  k_T = k_20*theta^(Tw(j,i)-273.15);  % F_t: temperature effect limitation
%  PAR_0 = 0.47*H_G_m(j);              % photosynthetically active radiation [w/m^2]
%  PAR_z = PAR_0.*exp(-k_e.*H(i));       % photosynthetically active radiation [w/m^2]
%  F_lp = PAR_z/(PAR_z+K_lp);    % F_l: phytoplankton light effect limitation
%  r_photo = r_max*k_T*F_lp;     % photosynthesis rate


% Oxygen loss via RESPIRATION
r_resp = (k_resp_max./H).*(c_x(i)./(KO2+c_x(i))).*c_x_bod(i)./(KBOD+c_x_bod(i));

% Oxygen exchange via GAS TRANSFER
k2 = sqrt(Dm.*q./(H.^3)).*1.0241.^((Tw(j,i)-273.15)-20);
r_gas = k2.*(csat(j,i) - c_x(i)) ;

% reaction terms for oxygen and BOD concentrations = Sum of all reaction processes
 dcdt_react_O2 = r_photo - r_resp + r_gas ;
 dcdt_react_bod = - r_resp ;

end
end