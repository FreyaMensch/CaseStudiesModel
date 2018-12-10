function dTdt_A = Advectionfun(n,q,dx,T,T_in)
%ADVECTIONFUN Summary of this function goes here
%   Detailed explanation goes here
 for  i=2:n
dTdt_A(i)=(q(i)./dx).*(T(i-1)-T(i));
  end
dTdt_A(1)=(q(1)./dx).*(T_in-T(1));

end

