function [ RHSV ] = RHS_V(i,j,dL,dH,Ren,V,FV,GV)
%Calculates and returns Right hand side of Eqs for given values
%For each call returns right hand side for single point

%V velocity
%RHSV=-(DDX(i,j,dL,FV)+DDY(i,j,dH,GV))+(DDXX(i,j,dL,V)+DDYY(i,j,dH,V))/Ren;
dFVdx=0.5*(FV(i,j+1)-FV(i,j-1))/dL; %dFV/dx
dGVdy=0.5*(GV(i+1,j)-GV(i-1,j))/dH;  %dUV/dy
d2Vdx2=(V(i,j-1)-2*V(i,j)+V(i,j+1))/dL^2;
d2Vdy2=(V(i-1,j)-2*V(i,j)+V(i+1,j))/dH^2;
RHSV=-(dFVdx+dGVdy)+(d2Vdx2+d2Vdy2)/Ren;  
end

