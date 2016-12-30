clc
clear
%In this code we try to solve unsteady heat transfer problem over
%rectangular domain Toward the steady Soltion.
%Right and Left Walls are kept at T0, T1 Top and bottom walls are insulated
%Left and  Right Boundries are kept in T0 and T1, respectively
%Top and Buttom Boundries are Zero Flux
%BCS Can be change in the Bcs Function
%By Mohammad Aghakhani, 2012

%  *************************Top******************************
%  *........................................................*
%  *........................................................* 
% Left--->T0......................................T1<---Right
%  *........................................................*
%  *........................................................*
%  ***********************Bottom*****************************

%-----------------Inputs-----------
L=2;   %Channel Lenght
H=1;   %Chnnel with
T0=300; %Left wall Temprature
T1=50;%Right wall temperature
t0=100;   %iniial Values
alpha=0.00023; %Thermal diffusivity
m=100;  % No. of points along Top & Bottom
n=100 ; %No. of point along Left & Right sides
IT=1;    %Current iteration No.
MIT=100000; %Maximum allowabe iteration
Dt=0.22; %time step
eps=1e-6; %error
errT=1000; %Error in two con. time step
%------------------------------------

%Coordinate of nodes
X=zeros(n,m);
Y=zeros(n,m);

%Solution varibles @ Current  iteration
T=zeros(n,m);

%Solution varibles @ Previous iteration
Told=zeros(n,m);

%Grid Genration
[X,Y,dL,dH]=Grid(m,n,L,H);
F0=alpha*Dt/(dL*dH);
fprintf(1,'Fourier Number  =  %2.5f\n',F0);
disp('Press any Key to continue');
pause
if(F0>=.25)
    disp('Fourier Number is too High, Please decrease Time Step');
    return;
end

%Initiate the solution
[T]=initiate(n,m,t0);

%Set boundry condition for first time
[T]=Bcs(n,m,T,T0,T1);
Told=T;
error=zeros(1,MIT);
%Begin Iteration
while((IT<MIT)&&(errT>eps))
    IT=IT+1;
    %shift solution from old iteration
    Told=T;
    [T]=Bcs(n,m,T,T0,T1);
        %Caluate Right hand side(aplha*D^T/DX^2+*D^T/DY^2)
        %RHST=alpha*(DDXX(n,m,dL,T)+DDYY(n,m,dH,T));
        %Explicit Euler Method
    for i=2:n-1
         for j=2:m-1
             DTDXX=(  T(i,j+1)-2.*T(i,j)+T(i,j-1)  )/(dL.*dL);
             DTDYY=(  T(i+1,j) -2*T(i,j)+T(i-1,j)  )/(dH.*dH);
             T(i,j)=T(i,j)+alpha *( DTDXX+DTDYY)*Dt;
         end
    end
       %Calculate Errors
    errT=max(max(abs((T-Told))));
    error(IT)=errT;
    fprintf(1,'IT=%i  Time=%2.6e  Error=%2.6e\n',IT,IT*Dt,errT);
end
if(errT<eps)
    fprintf(1,'Converged in %i Iterations',IT);
else
    disp('Maximum Iteration Number Reached');
    plot(1:IT,log10(error(1:IT)),'- r');
    xlabel('Iteration');
    ylabel('Log(error)');
    title('Convergence History');
    return;
end
%Postprocces
plot(1:IT,log10(error(1:IT)),'- r');
xlabel('Iteration');
ylabel('Log(error)');
title('Convergence History');


figure
[C1,h1] = contourf(T,30);
text_handle = clabel(C1,h1,'manual');
colorbar
xlabel('x')
ylabel('y')
title(strcat('Temperature Contour @',num2str(IT*Dt)));
drawnow

figure
surf(X,Y,T);
title(strcat('Temperature Surface @',num2str(IT*Dt)));
xlabel('x')
ylabel('y')
zlabel('Temperature')
axis fill;

