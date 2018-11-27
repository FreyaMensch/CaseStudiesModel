clear all
close all
clc
%% Example Script from Video with comments
L=100;     %wall thickness[m] => length of reach 
n =10;      %number of sumulation nodes  => nr of cells in each reach
T0=8;      %initial temperature
T1s=40;    %surface 1 temperature  => boundary condition (Dirichlet bc) @ beginning of river
T2s=20;    %surface 2 temperature  => boundary condition (Dirichlet bc) @ end of river

dx= L/n;            %      => distance between the cell centers

alpha=1;         %thermal diffusivity    => we have to replace this by parameters for our purpose

t_final=6000;         %simulation time [s]
dt=10;           %fixed time step [s]  
x=dx/2:dx:L-dx/2;   % it can be here instead of dx/2 write just dx => no this does the simulation at the center of each cell and that's what we want!?
                    % this is the spatial vector

T = ones(n,1)*T0;   % => creates matrix that will be input by calculated values of Temperature values
dTdt=zeros(n,1);    % => creates matrix that will be input by calculated values of time derivative
t=0:dt:t_final;     % => created the temporal vector
T_eq = ones(length(x),length(t));

for j=1:length(t)   % temporal iteration
    for i= 2:n-1    % spacial iteration
        dTdt(i)=alpha*(-(T(i)-T(i-1))/dx^2+(T(i+1)-T(i))/dx^2);
    end
    dTdt(1)=alpha*(-(T(1)-T1s)/dx^2+(T(2)-T(1))/dx^2);          % bc outside the spatial loop
    dTdt(n)=alpha*(-(T(n)-T(n-1))/dx^2+(T2s-T(n))/dx^2);        % bc outside the spatial loop
    T = T+dTdt*dt;
    T_eq(:,j) = T;
    
    %T_eq(j) = T(j,:);
    % collecting all values in a matrix:
    % t=0:dt(1):te;
    % T_eq(length(x),length(t))=T(x(i))
    
    figure(1)
    plot(x,T,'linewidth',3)
    axis([0 L 0 50])
    xlabel(' Distance(m)')
    ylabel('Temperature(\circC)')
    pause(0.1)
end




% %% Example Script from Video with comments
% L=0.1;     %wall thickness[m] => length of reach 
% n =10;      %number of sumulation nodes  => nr of cells in each reach
% T0=0;      %initial temperature
% T1s=40;    %surface 1 temperature  => boundary condition (Dirichlet bc) @ beginning of river
% T2s=20;    %surface 2 temperature  => boundary condition (Dirichlet bc) @ end of river
% 
% dx= L/n;            %      => distance between the cell centers
% 
% alpha=0.01;         %thermal diffusivity    => we have to replace this by parameters for our purpose
% 
% t_final=60;         %simulation time [s]
% dt=0.001;           %fixed time step [s]  
% x=dx/2:dx:L-dx/2;   % it can be here instead of dx/2 write just dx => no this does the simulation at the center of each cell and that's what we want!?
%                     % this is the spatial vector
% 
% T = ones(n,1)*T0;   % => creates matrix that will be input by calculated values of Temperature values
% dTdt=zeros(n,1);    % => creates matrix that will be input by calculated values of time derivative
% t=0:dt:t_final;     % => created the temporal vector
% 
% for j=1:length(t)   % temporal iteration
%     for i= 2:n-1    % spacial iteration
%         dTdt(i)=alpha*(-(T(i)-T(i-1))/dx^2+(T(i+1)-T(i))/dx^2);
%     end
%     dTdt(1)=alpha*(-(T(1)-T1s)/dx^2+(T(2)-T(1))/dx^2);
%     dTdt(n)=alpha*(-(T(n)-T(n-1))/dx^2+(T2s-T(n))/dx^2);
%     T = T+dTdt*dt;
%     figure(1)
%     plot(x,T,'linewidth',3)
%     axis([0 L 0 50])
%     xlabel(' Distance(m)')
%     ylabel('Temperature(\circC)')
%     pause(0.1)
% end