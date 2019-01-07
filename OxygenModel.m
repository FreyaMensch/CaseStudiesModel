%% This is a subscript of RiverQuality_MASTER
%% OxygenModel

%% ///////////// DATA IMPUT ////////////////
c_air = 20.*rand(length(t),1);    % [mg/L]=[g/m³] oxygen concentration in air over time(random) 

%% /////////// TIME LOOP TO CALCULATE OXYGEN CONCENTRATION ///////////////
%% ///////////////////(with Euler method) ///////////////////////////
c_in = 10;                 % bc: input concentration at top end of river
c_x = ones(1,n)*0;        % initial concentration in water (whole river)
c_xt = ones(length(t),n)*282; % concentration over time and space
dcdt_A = zeros(1,n); dcdt_D = zeros(1,n); dcdt_S = zeros(1,n);

for j=1:length(t)
%========================= ADVECTION ======================================
dcdt_A = Advectionfun(n,q,dx,c_x,c_in);
%========================= DISPERSION =====================================
dcdt_D = Dispersionfun(n,Dbulk,A_c,dx,c_x);
% %========================= EXTERNAL SOURCE ================================
% only template - NOT Adjusted yet!!
% for  i=pos:pos             % point source (inserted in only one cell)
%      W = rho_w*C_p*Q_s*T_s   ;
%     dcdt_S(i)=W./(rho_w*C_p*A_c(i)*dx);
% end
% %========================= REACTIONS & MASS TRANSFER ======================
% % only started - NOT working / finished yet!!
% dcdt_react=Reactionfun(t,n,H_G_m);

%========================= SUM OF ALL =====================================
c_x = c_x + dcdt_A*dt + dcdt_D*dt + dcdt_S*dt + dcdt_react*dt;
%T_eq = T_eq + dTdt_A*dt + dTdt_D*dt + dTdt_S*dt ;  % T_eq is calculated by setting all air-water fluxes = 0
% 
c_xt(j,:) = c_x;
%Tw_eq(j,:) = T_eq;
end
