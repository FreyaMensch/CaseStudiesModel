%% /////////////////// 1. RIVER GEOMETRY //////////////////////////////// 
%% COMMENTS on Geometry
% info: so far 3 reaches are computed, all have different hydraulic 
% parameters, this can be changed to more reaches if liked
% output parameters: (needed for further calculations) H, q, A_c, A_s with
% B_s (width @surface) for Dispersion

%% PARAMETERS to be adjusted
Q = 0.5;                              % volumetric flux [m³/s]       % eg Neckar in Tübingen ~ 30 m³/s
n_man = 0.03;               % Manning roughness coefficient
L_tot = 30000;                       % total length of river [m]
reach_nr = 3;                         % number of reaches [-]
L_reach = L_tot/reach_nr; 
S_river = [0.001 0.002 0.003];       % bottom slope [m/m]
%S_river = [0.001 0.002 0.004];         % bottom slope [m/m]
B_0 = [2 3 4];                       % bottom width [m]
%B_0 = [2.5 3 4];                       % bottom width [m] 
s_bank = [1 1 1];                    % slope of river banks [m/m] (dy/dy)
%s_bank = [1 2 0.5];                    % slope of river banks [m/m] (dy/dy)
H_0 = 0.1;                             % initial water depth in first reach

% idea just for fun: S_river= 0.001*rand(10,1); B_m=(1-0.1).*rand(8724,1)+0.1;
% Tried, and it works > ask freya for code
%%
%�������������������� SPATIAL DISCRETIZATION  �����������������������������
dx = 50;                % element length fix [m]
n = L_tot/dx ;        %  number of cells
x = dx/2:dx:L_tot-dx/2 ;  % central point of every cell
%��������������������������������������������������������������������������

% river geometry parameters in every cell
S_river = repmat(S_river,L_reach/dx, 1); S_river = S_river(:)';
B_0 = repmat(B_0,L_reach/dx,1); B_0 = B_0(:)';
s_bank = repmat(s_bank,L_reach/dx,1); s_bank = s_bank(:)';
          
%% CALCULATION of river geometry parameters
H = zeros(1,n);            % creates a vector for water depth [m]
A_c = zeros(1,n);          % creates a vector for cross sectional area [m²]
A_s = zeros(1,n);          % creates a vector for surface area (water-atmosphere interface) [m²]
q = zeros(1,n);            % creates a vector for specific discharge [m/s]

for i=1:n
    % Manning's equation solved for water depth H:    [QUAL2K maual eq(17)]
    if i == 1
        H(i)= (((Q*n_man)^(3/5))*(B_0(i)+2.*H_0*sqrt(s_bank(i).^2+1))^(2/5))/((S_river(i)^(3/10)).*(B_0(i)+s_bank(i).*H_0));
    else 
        H(i)= (((Q*n_man)^(3/5)).*(B_0(i)+2.*H(i)*sqrt(s_bank(i).^2+1)).^(2/5))/((S_river(i)^(3/10)).*(B_0(i)+s_bank(i).*H(i)));
    end    
A_c(i) = (B_0(i)+s_bank(i)*H(i))*H(i);  % cross-sectional area between elements [m²] 
b_s(i) = H(i)/s_bank(i);                % length of water surface over sloped bank [m]
B_s(i) = B_0(i) + 2*b_s(i);             % total width of river at water surface [m]
A_s(i) = B_s(i)*L_reach;                % surface area (water-atmosphere interface) [m²]
q(i) = Q./A_c(i);                       % specific discharge [m/s]
end


%% time it takes the water from start to end of river
%t_res = sum(L_reach./q);             % residence time of water parcel in river [seconds]
t_res = sum(L_reach./q)/3600;        % residence time of water parcel in river [hours]