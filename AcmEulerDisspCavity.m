clc
clear
%Calculates Steady Cavity Flow by Artifical comoressibilty Method
% And Adds Numerical Dissipation to Equations
%By Mohammad Aghakhani
%-----------------Input-----------
L=1;%Channel Lenght
H=1;%Channel with
Ren=10; %Reynonld Number
Beta=1; %Artifical comoressibilty coeffcient
m=30 ; % No. of points along channel walls
n=30 ; %No. of point along channel sides
MIT=1000000; %Maximum allowabe iteration
%Dt=.0005; %time step
CFL=0.25; %Courant Number
eps=1e-4; %error
err=zeros(1,MIT);
err(1)=1000; %Error in two con. time step
epsx=0.1; %Dissipation Coeficient in X Direction
epsy=0.1;  %Dissipation Coeficient in Y Direction 
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
[P,U,V ] = BcsCavity( n,m,P,U,V );
%Determine minimum allowabe time step size
Dt=CFL_Test(Beta,CFL,dL,dH,U,V);
%Begin Iteration
fprintf(1,'Maximum Allowable Time Step=%2.6e\n',Dt);
disp('Press any key')
pause
IT=1;
while((IT<MIT)&&(err(IT)>eps))
    IT=IT+1;
    %shift solution from old iteration
    Pold=P;
    Uold=U;
    Vold=V;
    %Apply bcs
    [P,U,V]=BcsCavity(n,m,P,U,V);

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
           %Add diissioation to RHS
            RHSP=RHSP-(  DisspX( i,j,m,P,epsx)+ DisspY( i,j,n,P,epsy)  );
            RHSU=RHSU-( DisspX( i,j,m,U,epsx)+ DisspY( i,j,n,U,epsy) );
            RHSV=RHSV-( DisspX( i,j,m,V,epsx)+ DisspY( i,j,n,V,epsy)  );
            %Calculate Solution in new time step
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
   fprintf(1,'IT=%i   Error=%2.6e\n',IT,err(IT));
end

if(err(IT-1)>1000)
    fprintf(1,'Iterations Diverged\n');
    fprintf(1,'Please Consider Change in Time Step Or other Parmeters and Run Code Again\n');
    disp('Press any key')   
    pause
end
 if(err(IT)<eps)
    fprintf(1,'Converged in %i Iterations\n',IT);
else
    disp('Maximum Iteration Number Reached');
    plot(1:IT,log10(err(1:IT)),'- r');
    xlabel('Iteration');
    ylabel('Log(error)');
    title('Convergence History');
    return;
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
text_handle = clabel(C1,h1);
colorbar
title('CONTOURES OF u(x,y)');
xlabel('x')
ylabel('y')
axis([0 L 0 H])
drawnow


figure
surf(X,Y,P)
title('Surface of Pressure')
xlabel('x')
ylabel('y')
zlabel('P')


figure
surf(X,Y,U)
title('Surface of U(x Velocity)')
xlabel('x')
ylabel('y')
zlabel('U')

figure
surf(X,Y,V)
title('Surface of V(y Velocity)')
xlabel('x')
ylabel('y')
zlabel('V')

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
hold off

figure
quiver(X,Y,U*10,V*10)
title('Velocity Vectors')
axis([0 L 0 H])

% figure
% streamline(X,Y,U,V,L/2,H/2);
% title('StreamLine Plot')
% axis([0 L 0 H])




