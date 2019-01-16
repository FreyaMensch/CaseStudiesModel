%% This is a subscript of RiverQuality_MASTER
%% OxygenModel

%% ///////////// DATA IMPUT ////////////////
load ('TemperatureModel_14_01_90days.mat','n','q','dx','x','t','dt','Dbulk','A_c','Dm','H','H_G_m','D')
load ('Tw_14_01_90days.mat')    %(Tw from temperature model)
c_air = 20.*rand(length(t),1);    % [mg/L]=[g/m³] oxygen concentration in air over time(random) 

%% /////////// TIME LOOP TO CALCULATE OXYGEN CONCENTRATION ///////////////
%% ///////////////////(with Euler method) ///////////////////////////
% saturation concentration [g/m3]
csat = 468./(31.6 + (Tw - 273.15));    %[g/m^3]

c_in = csat(1,1);                 % bc: input concentration at top end of river
c_x = zeros(1,n);       % initial concentration in water (whole river)
c_xt = ones(length(t),n); % concentration over time and space
dcdt_A = zeros(1,n); dcdt_D = zeros(1,n); dcdt_S = zeros(1,n);

c_in_bod = 10;                 % bc: input concentration at top end of river
c_x_bod = zeros(1,n);        % initial concentration in water (whole river)
c_xt_bod = ones(length(t),n); % concentration over time and space
dcdt_A_bod = zeros(1,n); dcdt_D = zeros(1,n); dcdt_S = zeros(1,n);
for j=1:length(t)
%========================= ADVECTION ======================================
dcdt_A = Advectionfun(n,q,dx,c_x,c_in);
dcdt_A_bod = Advectionfun(n,q,dx,c_x_bod,c_in_bod);
% ========================= DISPERSION =====================================
dcdt_D = Dispersionfun(n,Dbulk,A_c,dx,c_x);
dcdt_D_bod = Dispersionfun(n,Dbulk,A_c,dx,c_x_bod);

% %========================= EXTERNAL SOURCE ==============================
%  - NOT Adjusted yet!!
% for  i=pos:pos             % point source (inserted in only one cell)
%      W = Q_s*c_x   ;
%     dcdt_S(i)=W./(A_c(i)*dx);
% end
% %========================= REACTIONS & MASS TRANSFER ====================
% 
[dcdt_react_O2, dcdt_react_bod] = Reactionfun(t,j,n,H_G_m,c_x,c_x_bod,H,Dm,q,Tw,csat);

%========================= SUM OF ALL =====================================
c_x = c_x + dcdt_A*dt + dcdt_D*dt + dcdt_react_O2*dt ;
c_x_bod = c_x_bod + dcdt_A_bod*dt + dcdt_D_bod*dt + dcdt_react_bod*dt ;

% + dcdt_S*dt

% figure (1)
% plot(x,c_x),hold on
% plot(x,c_x_bod)
% ylim ([0 20])
% xlim ([0 30000])
% hold off
c_xt(j,:) = c_x;
c_xt_bod(j,:) = c_x_bod ;
end

%%
% plot(x,c_xt(120000,:)), hold on;
% plot(x,c_xt_bod(120000,:))
% legend