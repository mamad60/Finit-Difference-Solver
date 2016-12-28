function [ RHSP ] = RHS_P(i,j,dL,dH,FP,GP)
%Calculates and returns Right hand side of Eqs for given values
%For each call returns right hand side for single point

%Pressure
%RHSP=-(DDX(i,j,dL,FP)+DDY(i,j,dH,GP))+0+0;
dFPdx=0.5*(FP(i,j+1)-FP(i,j-1))/dL;   %dFP/dx
dGPdy=0.5*(GP(i+1,j)-GP(i-1,j))/dH;   %dGP/dy
RHSP=-(dFPdx+dGPdy);
end


