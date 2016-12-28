clc
clear
%Solves Channel FLow by Artifical comoressibilty coeffcient 
%and 4th Order Rung-Kutta Method
%And Adds Numerical Dissipation
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
   
    %Shift solution from old iteration U(0) of Rung-kutta
    Pold=P;
    Uold=U;
    Vold=V; %U(0)
    %Apply Bcs
    [P,U,V]=Bcs(n,m,P,U,V);
 [FP,GP,FU,GU,FV,GV] = Flux(P,U,V,Beta);

        
    %4th Oreder Rung-Kutta Method
    for i=2:n-1
         for j=2:m-1
%------------------Calculation of A coeffients--------
              %Compute Right hand sides of Eqs PU(0)
            [ RHSUA ] = RHS_U(i,j,dL,dH,Ren,U,FU,GU);
            [ RHSVA ] = RHS_V(i,j,dL,dH,Ren,V,FV,GV);
            [ RHSPA ] = RHS_P(i,j,dL,dH,FP,GP);
            %Add diissioation to RHS  P(0)
            RHSUA=RHSUA-( DisspX( i,j,m,U)+ DisspY( i,j,n,U) );
            RHSVA=RHSVA-( DisspX( i,j,m,V)+ DisspY( i,j,n,V)  );
            RHSPA=RHSPA-(  DisspX( i,j,m,P)+ DisspY( i,j,n,P)  );
            %Calculate Solution in new time step U(1)
            U(i,j)=Uold(i,j)+RHSUA*Dt/2;
            V(i,j)=Vold(i,j)+RHSVA*Dt/2;
            P(i,j)=Pold(i,j)+RHSPA*Dt/2;
           
 %------------------Calculation of B coeffients--------
            [FP(i,j),GP(i,j),FU(i,j),GU(i,j),FV(i,j),GV(i,j)] = Fluxij(i,j,P,U,V,Beta);     
              %Compute Right hand sides of Eqs
            [ RHSUB ] = RHS_U(i,j,dL,dH,Ren,U,FU,GU);
            [ RHSVB ] = RHS_V(i,j,dL,dH,Ren,V,FV,GV);
            [ RHSPB ] = RHS_P(i,j,dL,dH,FP,GP);  % PU(1)
            %Add diissioation to RHS
            RHSUB=RHSUB-( DisspX( i,j,m,U)+ DisspY( i,j,n,U) );
            RHSVB=RHSVB-( DisspX( i,j,m,V)+ DisspY( i,j,n,V)  );
            RHSPB=RHSPB-(  DisspX( i,j,m,P)+ DisspY( i,j,n,P)  );
            %Calculate Solution in new time step
            U(i,j)=Uold(i,j)+RHSUB*Dt/2;
            V(i,j)=Vold(i,j)+RHSVB*Dt/2;
            P(i,j)=Pold(i,j)+RHSPB*Dt/2;  %U(2)
 %------------------Calculation of C   coeffients--------
            [FP(i,j),GP(i,j),FU(i,j),GU(i,j),FV(i,j),GV(i,j)] = Fluxij(i,j,P,U,V,Beta);     
              %Compute Right hand sides of Eqs
            [ RHSUC ] = RHS_U(i,j,dL,dH,Ren,U,FU,GU);
            [ RHSVC ] = RHS_V(i,j,dL,dH,Ren,V,FV,GV);
            [ RHSPC ] = RHS_P(i,j,dL,dH,FP,GP);  % PU(2)
            %Add diissioation to RHS
            RHSUC=RHSUC-( DisspX( i,j,m,U)+ DisspY( i,j,n,U) );
            RHSVC=RHSVC-( DisspX( i,j,m,V)+ DisspY( i,j,n,V)  );
            RHSPC=RHSPC-(  DisspX( i,j,m,P)+ DisspY( i,j,n,P)  );
            %Calculate Solution in new time step
            U(i,j)=Uold(i,j)+RHSUC*Dt;
            V(i,j)=Vold(i,j)+RHSVC*Dt;
            P(i,j)=Pold(i,j)+RHSPC*Dt;  %U(3)
%------------------Calculation of D coeffients--------
            [FP(i,j),GP(i,j),FU(i,j),GU(i,j),FV(i,j),GV(i,j)] = Fluxij(i,j,P,U,V,Beta);     
              %Compute Right hand sides of Eqs
            [ RHSPD ] = RHS_P(i,j,dL,dH,FP,GP);  % PU(3)
            [ RHSUD ] = RHS_U(i,j,dL,dH,Ren,U,FU,GU);
            [ RHSVD ] = RHS_V(i,j,dL,dH,Ren,V,FV,GV); 
            %Add diissioation to RHS
            RHSUD=RHSUD-( DisspX( i,j,m,U)+ DisspY( i,j,n,U) );
            RHSVD=RHSVD-( DisspX( i,j,m,V)+ DisspY( i,j,n,V)  );
            RHSPD=RHSPD-(  DisspX( i,j,m,P)+ DisspY( i,j,n,P)  );

%                %U(4) New Time Step Calculations
            U(i,j)=Uold(i,j)+(RHSUA+2*RHSUB+2*RHSUC+RHSUD)*Dt/6;   
            P(i,j)=Pold(i,j)+(RHSVA+2*RHSVB+2*RHSVC+RHSVD)*Dt/6;   
            P(i,j)=Pold(i,j)+(RHSPA+2*RHSPB+2*RHSPC+RHSPD)*Dt/6;   
            

            
            
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

