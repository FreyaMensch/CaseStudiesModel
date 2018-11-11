% function [Dr,Dl]=Dispersion(time, q,B,Hmean,So,Ac)
 Vmean= 0.1;   %  mean velocity [m/s]  from Freya
 w=4   ;    %  width of river [m]   Which one???? from Freya
 Hmean=1 ;   % mean depth [m]  from Freya
 So= 0.01 ;  %    channel slope [-]   data input in Freya
 g=9.81 ;   %[m/s^2]
Ac=4  ;        %cross-sectional area   from Freya
% dx=100;        distance between elements [m]
L=10000;  % reach length [m]
esp=0.01;    %delta x [m]
x=[0:esp:L];
dx=x(2:end)-x(1:end-1);
Vsh=sqrt(g*dmean^2*So);    %% Shear velocity [m/s](QUAL2K.p18) depth will change in every reach?? 

Dnum=Vmean*dx/2 ;  % numerical dispersion right side (QUAL2K.p18)
 
D=0.011*Vmean^2*w^2/(dmean*Vsh)  ;%%  (QUAL2K.p18) check which variables are going to change in every reach
 
% if Dnum<=D
%     Dm=D-Dnum;
% else if Dnum>D
%     Dm=0;
%     end
% end
% 
% % CONSIDER BOUNDARY CONDITIONS  (ZERO DISPERSION AND DIRICHLET)!!!!!
% 
 Dbulk= D*Ac*2./(dx(1:end-1)+dx(2:end))  ;
 
%  Dbulk= Dm*Ac*2./(dx(1:end-1)+dx(2:end))  ;
 
 

% end
 