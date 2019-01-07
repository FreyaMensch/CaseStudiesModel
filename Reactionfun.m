function dcdt_react = Reactionfun(t,n,H_G_m)
%REACTIONFUN Summary of this function goes here
%   Detailed explanation goes here

for  i=1:n
% Oxygen input via PHOTOSYNTHESIS
    k_photo = 10;       % has to be made Temperature dependant!!
    r_photo = k_photo*H_G_m;

% Oxygen loss via RESPIRATION
    r_resp = 

% Oxygen exchange via GAS TRANSFER
    r_gasTrans =
    
    
% total change of oxygen concentration  over time = Sum of all processes
 dcdt_react(i) = 1;
 %dTdt_R(i) = 1./(H(i)*rho_w*C_p).*(H_insw(i)- H_outsw(i)- H_lwout(i)- H_ET(i) - H_c(i)+ H_lwin(i));
end
end