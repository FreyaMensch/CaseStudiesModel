function [Dbulk] = BulkDispersion(q,H,S_river,w,dx,Ac)
% Calculates bulk longitudinal dispersion. Input arguments:
% velocity q% [m/s], Depths of river [m], Bottom river slope [-], width of
% river [m], element lentgh [m], Cross-sectional area [m^2]
Vmean= q;   %  mean velocity [m/s]  from above
Hmean=H(1:3) ;     % mean depth [m]  from Freya
g=9.81 ;      %[m/s^2]
Vsh=sqrt((g*Hmean.^2).*S_river);    %% Shear velocity [m/s](QUAL2K.p18)
D=(0.011.*(Vmean.^2).*(w.^2))./(Hmean.*Vsh) ;%%  (QUAL2K.p18)

% dx=round(2*D/q,-1)    % element length ; Neumman number = 1/3 , Courant number = 1!!!

% dx=x(2:end)-x(1:end-1);
Dnum=Vmean*dx/2 ;  % numerical dispersion (QUAL2K.p18) 
if Dnum<=D
    Dm=D-Dnum;
else if Dnum>D
    Dm=0;
    end
end

Dbulk= (Dm.*Ac.*2)./dx  ;   % assuming that dx does not change
end

