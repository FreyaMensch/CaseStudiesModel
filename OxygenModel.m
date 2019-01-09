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
% %========================= EXTERNAL SOURCE ==============================
% only template - NOT Adjusted yet!!
% for  i=pos:pos             % point source (inserted in only one cell)
%      W = rho_w*C_p*Q_s*T_s   ;
%     dcdt_S(i)=W./(rho_w*C_p*A_c(i)*dx);
% end
% %========================= REACTIONS & MASS TRANSFER ====================
% % only started - NOT working / finished yet!!
% dcdt_react=Reactionfun(t,n,H_G_m);

    % photosynthesis rate
    %F_t=k_20*theta^(T-20);       % Ft: temperature effect limitation
    
   
    theta=1.056 % in general
    k_20=       % constant at 20 Celsius
    ke=         % light extintion coeff [1/m]
    H=0.1       % water depth
    K_lb=       % bottom algae parameter [w/m2]
    PAR=PAR_0*(exp(1)^(-ke*H));   % photosynthetically active radiation
    PAR_0=0.47*H_G;  % photosynthetically active radiation [w/m^2]
    
    
    F_t=k_20*theta^(T-273.15);        % F_t: temperature effect limitation
    F_l=PAR/(PAR+K_lb);               % F_l: light effect limitation

    r_photo=r_max*F_t*F_l;            % photosynthesis rate


    % respiration
    %r_resp=k_resp_max/h.*c(:,1)./(KO2+c(:,1)).*c(:,2)./(KBOD+c(:,2));%
    %why divided by h?
    r_resp=k_resp_max.*(c(:,1)./(KO2+c(:,1)).*c(:,2)./(KBOD+c(:,2))); % need to consider consenctration
   
    % gas transfer
    k2 = 1/h*sqrt(Dm*v/h)*1.0241.^(T-20);
    %k2 = sqrt(Dm*v/(h^3))*1.0241.^(T-20);



%========================= SUM OF ALL =====================================
c_x = c_x + dcdt_A*dt + dcdt_D*dt + dcdt_S*dt + dcdt_react*dt;
%T_eq = T_eq + dTdt_A*dt + dTdt_D*dt + dTdt_S*dt ;  % T_eq is calculated by setting all air-water fluxes = 0
% 
c_xt(j,:) = c_x;
%Tw_eq(j,:) = T_eq;
end
