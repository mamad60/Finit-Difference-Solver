function [ RHSP,RHSU,RHSV ] = RHS(i,j,dL,dH,Ren,U,V,FP,FU,FV,GP,GU,GV)
%Calculates and returns Right hand side of Eqs for given values
%For each call returns right hand side for single point

%Pressure
%RHSP=-(DDX(i,j,dL,FP)+DDY(i,j,dH,GP))+0+0;
dFPdx=0.5*(FP(i,j+1)-FP(i,j-1))/dL;   %dFP/dx
dGPdy=0.5*(GP(i+1,j)-GP(i-1,j))/dH;   %dGP/dy
RHSP=-(dFPdx+dGPdy);

%U velocity
%RHSU=-(DDX(i,j,dL,FU)+DDY(i,j,dH,GU))+(DDXX(i,j,dL,U)+DDYY(i,j,dH,U))/Ren;
dFUdx=0.5*(FU(i,j+1)-FU(i,j-1))/dL; %dFU/dx
dGUdy=0.5*(GU(i+1,j)-GU(i-1,j))/dH;  %dGU/dy
d2Udx2=(U(i,j-1)-2*U(i,j)+U(i,j+1))/dL^2;%
d2Udy2=(U(i-1,j)-2*U(i,j)+U(i+1,j))/dH^2;
RHSU=-(dFUdx+dGUdy)+(d2Udx2+d2Udy2)/Ren;

%V velocity
%RHSV=-(DDX(i,j,dL,FV)+DDY(i,j,dH,GV))+(DDXX(i,j,dL,V)+DDYY(i,j,dH,V))/Ren;
dFVdx=0.5*(FV(i,j+1)-FV(i,j-1))/dL; %dFV/dx
dGVdy=0.5*(GV(i+1,j)-GV(i-1,j))/dH;  %dUV/dy
d2Vdx2=(V(i,j-1)-2*V(i,j)+V(i,j+1))/dL^2;
d2Vdy2=(V(i-1,j)-2*V(i,j)+V(i+1,j))/dH^2;
RHSV=-(dFVdx+dGVdy)+(d2Vdx2+d2Vdy2)/Ren;  
end

