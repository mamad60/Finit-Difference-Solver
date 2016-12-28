function [FP,GP,FU,GU,FV,GV] = Fluxij(i,j,P,U,V,Beta)
%This function returns fluxes in Eq. as needed

%pressure    
FP=(Beta^2)*U(i,j);
GP=(Beta^2)*V(i,j);
%U
FU=U(i,j)^2+P(i,j);
GU=U(i,j)*V(i,j);
%V
FV=U(i,j)*V(i,j);
GV=V(i,j)^2+P(i,j);
 
end

