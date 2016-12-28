function [ Dt ] = CFL_Test(Beta,CFL,dL,dH,U,V )
%Calculates Maxium allowable Time step;

C=sqrt(U.^2+V.^2)+sqrt(Beta); 
Cmax=max(max(C)); %maximum sound+fluid speed
    
Lmin=min(dH,dL); %minimum grid spacing
Dt=CFL*Lmin/Cmax;
end

