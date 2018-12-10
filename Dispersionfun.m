function dTdt_D = Dispersionfun(n,Dbulk,A_c,dx,T)
%DISPERSIONFUN Summary of this function goes here
%   Detailed explanation goes here
 for  i=2:n-1
dTdt_D(i)=(Dbulk(i-1)./(A_c(i)*dx)).*(T(i-1)-T(i)) + (Dbulk(i)./(A_c(i)*dx)).*(T(i+1)-T(i)) ;
  end
dTdt_D(1)=0 ; % boundary conditions
dTdt_D(n)=0 ;
end

