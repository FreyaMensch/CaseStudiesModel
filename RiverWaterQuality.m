% Case Studies 2018-19
% River Water Quality Model
% by Veronica, Alireza and Freya
% last updated on 27.10.18 at 19:00 by Freya

%% To avoid conflicts: Who is working on which part?
% River geometry: Freya
% Radiation Balance: Alireza
% Oxygen model: Veronica


%%
%% MODEL
%%


%% ---------- RIVER GEOMETRY ----------- 
%% COMMENTS on Geometry
% info: so far 3 reaches are computed, all have different hydraulic 
% parameters, this can be changed to more reaches if liked
% output parameters: (needed for further calculations) H, A_c, A_s, q

%% PARAMETERS to be adjusted
Q = 10;                              % volumetric flux [m³/s]       % eg Neckar in Tübingen ~ 30 m³/s
n = 0.03;                            % Manning roughness coefficient
 S_river = [0.001 0.001 0.001];       % bottom slope [m/m]
%S_river = [0.001 0.002 0.004];       % bottom slope [m/m]
 B_0 = [3 3 3];                       % bottom width [m]
%B_0 = [2.5 3 4]; 
s_bank = [1 1 1];                    % slope of river banks [m/m] (dy/dy)
% s_bank = [1 2 0.5];
H_0 = 2;                             % initial water depth in first reach


%% CALCULATION of water depth
H = [H_0; zeros(3,1)]';     % creates a vector for water depth [m]
A_c = zeros(3,1)';          % creates a vector for cross sectional area [m²]
A_s = zeros(3,1)';          % creates a vector for surface area (water-atmosphere interface) [m²]
q = zeros(3,1)';            % creates a vector for specific discharge [m/s]

for j=1:3
for i=2:4
    % Manning's equation solved for water depth H:    [QUAL2K maual eq(17)]
    %H(i)= (((Q*n)^(3/5))*(B_0(i-1)+2.*H(i-1)*sqrt(s_bank(i-1).^2+1))^(2/5))/((S_river(i-1)^(3/10)).*(B_0(i-1)+s_bank(i-1).*H(i-1)));
    % same as above, BUT all H refering to initial height H_0:
    H(i)= (((Q*n)^(3/5))*(B_0(i-1)+2.*H_0*sqrt(s_bank(i-1).^2+1))^(2/5))/((S_river(i-1)^(3/10)).*(B_0(i-1)+s_bank(i-1).*H_0));
% Note: the first iteration of H(i-1) refers to the H_0, the following to 
% the H value calculated in the previous iteration whereas the other 
% parameters(i-1) all refer to the parameters defined in PARAMETERS section
end
A_c(j) = (B_0(j)+s_bank(j)*H(j+1))*H(j+1);     % cross-sectional area [m²] 
%A_s(j) = ;                                     % surface area (water-atmosphere interface) [m²]
q(j) = Q./A_c(j);                              % specific discharge [m/s]
end


%%
%% Longitudinal Dispersion      Manual p.18
% parameters
U

% longitudinal dispersion
E_p = 0.011*((U_i^2*B_i^2)/(H_i^2*Us_i^2));         % Fischer et al 1979; Manual p.18




%% HEAT SURFACE FLUX
%% Solar Radiation




%% Atmospheric Long-wave radiation

%%% change something


%% Water-longwave radiation



%% Conduction and Convection



























