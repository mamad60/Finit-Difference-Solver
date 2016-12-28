clc
clear
%This Calculates of Channel Flow with Artifical comoressibilty coeffcient
%And Adds 4th order Dissipation to Equations
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
eps=1e-4; %error
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
    
    %Explicit Euler Method
    for i=2:n-1
         for j=2:m-1
            %Compute Right hand sides of Eqs
            [ RHSP ] = RHS_P(i,j,dL,dH,FP,GP);
            [ RHSU ] = RHS_U(i,j,dL,dH,Ren,U,FU,GU);
            [ RHSV ] = RHS_V(i,j,dL,dH,Ren,V,FV,GV);
            %Add diissioation to RHS
            RHSP=RHSP-(  DisspX( i,j,m,P)+ DisspY( i,j,n,P)  );
            RHSU=RHSU-( DisspX( i,j,m,U)+ DisspY( i,j,n,U) );
            RHSV=RHSV-( DisspX( i,j,m,V)+ DisspY( i,j,n,V)  );
            %Calculate Solution in new time step Using Explicit Euler
            [P(i,j)]=Euler(P(i,j),Dt,RHSP);
            [U(i,j)]=Euler(U(i,j),Dt,RHSU);
            [V(i,j)]=Euler(V(i,j),Dt,RHSV);
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

