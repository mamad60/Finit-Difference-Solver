function [ RHSU ] = RHS_U(i,j,dL,dH,Ren,U,FU,GU)
%Calculates and returns Right hand side of Eqs for given values
%For each call returns right hand side for single point
%U velocity
%RHSU=-(DDX(i,j,dL,FU)+DDY(i,j,dH,GU))+(DDXX(i,j,dL,U)+DDYY(i,j,dH,U))/Ren;
dFUdx=0.5*(FU(i,j+1)-FU(i,j-1))/dL; %dFU/dx
dGUdy=0.5*(GU(i+1,j)-GU(i-1,j))/dH;  %dGU/dy
d2Udx2=(U(i,j-1)-2*U(i,j)+U(i,j+1))/dL^2;%
d2Udy2=(U(i-1,j)-2*U(i,j)+U(i+1,j))/dH^2;
RHSU=-(dFUdx+dGUdy)+(d2Udx2+d2Udy2)/Ren;
end

