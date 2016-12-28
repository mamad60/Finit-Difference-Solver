function [FP,GP,FU,GU,FV,GV] = Flux(P,U,V,Beta)
%This function returns fluxes in Eq. as needed

%pressure    
FP=(Beta^2)*U;
GP=(Beta^2)*V;
%U
FU=U.^2+P;
GU=U.*V;
%V
FV=U.*V;
GV=V.^2+P;
 
end

