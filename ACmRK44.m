clc
clear

%-----------------Input-----------
L=8;%Channel Lenght
H=1;%Chnnel with
Ren=100; %Reynonld Number
Beta=1.2; %Artifical comoressibilty coeffcient
m=60 ; % No. of points along channel walls
n=30 ; %No. of point along channel sides
MIT=100000; %Maximum allowabe iteration
 Dt=.005; %time step
%CFL=0.1; %Courant Number
eps=1e-6; %error
err=zeros(1,MIT);
err(1)=1000; %Error in two con. time step
%------------------------------------

%Coordinate of nodes
X=zeros(n,m);
Y=zeros(n,m);

%Solution varibles @ Current  iteration
U=zeros(n,m);
V=zeros(n,m);
P=zeros(n,m);
%intermediate varables
Um=zeros(n,m);
Vm=zeros(n,m);
Pm=zeros(n,m);

%Solution varibles @ Previous iteration
Uold=zeros(n,m);
Vold=zeros(n,m);
Pold=zeros(n,m);


%Flux & Resiuduals
    %Preesure
FP=zeros(n,m);
GP=zeros(n,m);

    %U
FU=zeros(n,m);
GU=zeros(n,m);
    %V
FV=zeros(n,m);
GV=zeros(n,m);
%Right hand sites of Eqs
RHSP=0;
RHSU=0;
RHSV=0;


%Intermediate Flux & Resiuduals
    %Preesurem
FPm=zeros(n,m);
GPm=zeros(n,m);

    %Um
FUm=zeros(n,m);
GUm=zeros(n,m);
    %Vm
FVm=zeros(n,m);
GVm=zeros(n,m);



%Grid Genration
[X,Y,dL,dH]=Grid(m,n,L,H);


%iniial Guess
u0=1;
v0=0;
p0=0.001;
%Initiate the solution
[P,U,V]=initiate(n,m,p0,u0,v0);
[P,U,V ] = Bcs( n,m,P,U,V );
%Determine minimum allowabe time step size
%Dt=CFL_Test(Beta,CFL,dL,dH,U,V);
%Begin Iteration

%disp('----------------------------------------------------------')
%disp('   IT        err(IT)      errP        errU       errV')
IT=1;
while((IT<=MIT)&&(err(IT)>eps))
   
    %shift solution from old iteration
    Pold=P;
    Uold=U;
    Vold=V;
    %Apply bcs
    [P,U,V]=Bcs(n,m,P,U,V);

    %Compute Fluxes From initial values From P,V,T
    [FP,GP,FU,GU,FV,GV] = Flux(P,U,V,Beta);
  %-------------------Set Intermediate Fluxes For used by -RK4-- 
    FPm=FP;
    GPm=GP;
    FUm=FU;
    GUm=GU;
    FVm=FV;
    GVm=GV;
    %-------------------------------------------------------------------------------
    
    %4th Order Rung Kutta
    for i=2:n-1
         for j=2:m-1
  %************************************************************************

 %----------------------------k1--------------------------------------------
           
            %----------U-------------------

            [ RHSU ] = RHS_U(i,j,dL,dH,Ren,U,FU,GU);
            k1u=Dt*RHSU;
            %----------V-------------------

            [ RHSV ] = RHS_V(i,j,dL,dH,Ren,V,FV,GV);
            k1v=Dt*RHSV;
            %----------P-------------------
 
            [ RHSP ] = RHS_P(i,j,dL,dH,FP,GP);
            k1p=Dt*RHSP;
            
 %----------------------------k2--------------------------------------------
            
             %----------U-------------------
             
            Um(i,j)=U(i,j)+k1u/2;
            FUm(i,j)= F_U(P(i,j),Um(i,j));
            GUm(i,j)=G_U(Um(i,j),V(i,j));
                  
            [ RHSU ] = RHS_U(i,j,dL,dH,Ren,Um,FUm,GUm);
            k2u=Dt*RHSU;
            
             %----------V-------------------
            
            Vm(i,j)=V(i,j)+k1v/2;
            FVm(i,j)=F_V(U(i,j),Vm(i,j));
            GVm(i,j)=G_V(P(i,j),Vm(i,j));
                      
            [ RHSV ] = RHS_V(i,j,dL,dH,Ren,Vm,FVm,GVm);
            k2v=Dt*RHSV;
            
             %----------P------------------
                      
             
            FPm(i,j)= F_P(Um(i,j),Beta);
            GPm(i,j)=G_P(Vm(i,j),Beta);
            
            [ RHSP ] = RHS_P(i,j,dL,dH,FPm,GPm);
            k2p=Dt*RHSP;
            
%----------------------------k3--------------------------------------------
            
            %----------U-------------------
            
            Um(i,j)=U(i,j)+k2u/2;
            FUm(i,j)= F_U(P(i,j),Um(i,j));
            GUm(i,j)=G_U(Um(i,j),V(i,j));                     
                      
            [ RHSU ] = RHS_U(i,j,dL,dH,Ren,Um,FUm,GUm);
            k3u=Dt*RHSU;
            
             %----------V-------------------
            
            Vm(i,j)=V(i,j)+k2v/2;         
            FVm(i,j)=F_V(U(i,j),Vm(i,j));
            GVm(i,j)=G_V(P(i,j),Vm(i,j));
            
            [ RHSV ] = RHS_V(i,j,dL,dH,Ren,Vm,FVm,GVm);
            k3v=Dt*RHSV;
            
             %----------P------------------
             
            FPm(i,j)= F_P(Um(i,j),Beta);
            GPm(i,j)=G_P(Vm(i,j),Beta);
            
            [ RHSP ] = RHS_P(i,j,dL,dH,FPm,GPm);
            k3p=Dt*RHSP;

%----------------------------k4--------------------------------------------
            %----------U-------------------
            
            Um(i,j)=U(i,j)+k3u; 
            FUm(i,j)= F_U(P(i,j),Um(i,j));
            GUm(i,j)=G_U(Um(i,j),V(i,j));
              
            [ RHSU ] = RHS_U(i,j,dL,dH,Ren,Um,FUm,GUm);
            k4u=Dt*RHSU;
            
             %----------V-------------------
             
            Vm(i,j)=V(i,j)+k3v;         
            FVm(i,j)=F_V(U(i,j),Vm(i,j));
            GVm(i,j)=G_V(P(i,j),Vm(i,j));
                        
            [ RHSV ] = RHS_V(i,j,dL,dH,Ren,Vm,FVm,GVm);
            k4v=Dt*RHSV;
              %----------P------------------
             
            FPm(i,j)= F_P(Um(i,j),Beta);
            GPm(i,j)=G_P(Vm(i,j),Beta);
            
            [ RHSP ] = RHS_P(i,j,dL,dH,FPm,GPm);
            k4p=Dt*RHSP;        
                   
                     
 %--------------------------Calculate Solution in New time Step---------------------------------------------

             
             
            U(i,j)=U(i,j)+(k1u+2*k2u+2*k3u+k4u)/6;
            V(i,j)=V(i,j)+(k1v+2*k2v+2*k3v+k4v)/6;  
            P(i,j)=P(i,j)+(k1p+2*k2p+2*k3p+k4p)/6; 
             
             
%****************************************************************************************            

         end
    end
     
   %Calculate Error
    errP=max(max(abs((P-Pold))));
    errU=max(max(abs((U-Uold))));
    errV=max(max(abs((V-Vold))));
    erri=max(errU,errV);
    err(IT)=max(erri,errP);
   
 IT=IT+1;
 err(IT)=err(IT-1);
 %[IT err(IT)  errP  errU  errV]
 %[IT err(IT)]
 fprintf(1,'IT=%i   Error=%2.6e\n',IT,err(IT));
end
disp('****************************************************************')
fprintf(1,'Beta=%2.2f       IT=%u       Dt=%2.2e        IT*Dt=%2.2f    Beta*Dt=%2.5f\n ',Beta, IT, Dt, IT*Dt, Beta*Dt);
disp('Press any key')
pause
Min=sum((U(:,1))*(dH));
Mout=sum((U(:,m-1))*(dH));
Mass_inbalance=Mout-Min
figure
plot( 1:IT,log10(err(1:IT)),'-. g');
xlabel('Iteration')
ylabel('Log10(Error)')
title('Error History')

figure
[C1,h1] = contourf(U,20);
text_handle = clabel(C1,h1,'manual');
colorbar
title('CONTOURES OF u(x,y)');
xlabel('x(1:n)')
ylabel('y(1:m)')

figure
surf(X,Y,P)
title('Surface of Pressure')

figure
surf(X,Y,U)
title('Surface of U(x Velocity)')


figure
surf(X,Y,V)
title('Surface of V(y Velocity)')


figure
hold on
plot(U(1:n,1),Y(1:n,1),'-* r')
plot(U(1:n,m/10),Y(1:n,m/10),'- g')
plot(U(1:n,3*m/10),Y(1:n,3*m/10),'-S m')
plot(U(1:n,5*m/10),Y(1:n,5*m/10),'- y')
plot(U(1:n,6*m/10),Y(1:n,6*m/10),'- b')
plot(U(1:n,7*m/10),Y(1:n,7*m/10),'-. c')
plot(U(1:n,8*m/10),Y(1:n,8*m/10),'- k')
plot(U(1:n,m),Y(1:n,m),'-- r')
legend('x=0','x=m/10','x=3*m/10','x=5*m/10','x=6*m/10','x=7*m/10','x=8*m/10','x=m',1)
title('Profiles of the U velocity across channel Section')

figure
quiver(X,Y,U*3,V*3)
time=IT*Dt
title('Velocity Vectors')

