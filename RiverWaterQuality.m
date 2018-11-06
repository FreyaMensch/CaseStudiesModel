% Case Studies 2018-19
% River Water Quality Model
% by Veronica, Alireza and Freya
% last updated on 27.10.18 at 19:00 by Freya


% TEST 2 and 3
% last try
% test ali
%test2 ali

%% FLOW BALANCE
%% Flow Balance
Q =                 % total flow rate [m³/s] - if everywhere the same
% Neckar in Tübingen ~ 30 m³/s

%% Mannings equation, River geometry
% parameters
S_river = 0.002;       % bottom slope [m/m]
n = 0.03;              % Manning roughness coefficient
B_0 = 3;               % bottom width [m]
s_bank = 1;          % slope of river banks [m/m]
H = 2;                  % water depth [m]

% cross-sectional area [m²]
A_c = (B_0+s_bank*H)*H; 

% wetted perimeter [m]
P = B_0 + 2*H*sqrt(s_bank^2+1);

% Mannings equation
Q = (sqrt(S_river)*A_c^(5/3))/(n*P^(2/3))          % total flow rate [m³/s]



%% Longitudinal Dispersion      Manual p.18
% parameters
U

% longitudinal dispersion
E_p = 0.011*((U_i^2*B_i^2)/(H_i^2*Us_i^2));         % Fischer et al 1979; Manual p.18




%% HEAT SURFACE FLUX
%% Solar Radiation

% Test Freya2 xxxxxxxx
% Test 3 freya

%% Atmospheric Long-wave radiation




%% Water-longwave radiation



%% Conduction and Convection



%% Sediment-Water Heat transfer























