function dTdt = AdvDisFun(~,T,x,Q,V,Dbulk)
%ADVFUN3 Summary of this function goes here
%   Detailed explanation goes here

%=========================== ADVECTION ===================================

dTdt_A = zeros(length(V),1);
dTdt_A(2:end) = (Q/V(1)).*(T(1:end-1)-T(2:end));
dTdt_A(1) = dTdt_A(1);

%=========================== DISPERSION ===================================

dTdt_D = zeros(length(V),1);
Hd = Dbulk(1)/V(1).*(T(1:end-1)-T(2:end));
dTdt_D = Hd(1:end-1)-Hd(2:end);
dTdt_D = [0; dTdt_D ;0];


%=========================== RADIATION ===================================

%  dTdt = dTdt_A;
dTdt = dTdt_A + dTdt_D;
end

