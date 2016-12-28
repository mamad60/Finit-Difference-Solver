clc
clear
%In this code we try to solve unsteady heat transfer problem over
%rectangular domain
%Right and Left Walls are kept at T0, T1 Top and bottom walls are insulated
%Left and  Right Boundries are kept in T0 and T1, respectively
%Top and Buttom Boundries are Zero Flux
%BCS Can be change in the Bcs Function

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
IT=0;    %Current iteration No.
MIT=5000; %Maximum allowabe iteration
Dt=0.2; %time step
eps=1e-5; %error
err=1000; %Error in two con. time step
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
if(F0>=.25)
    fprintf(1,'F0= %2.2f \n',F0)
end

%Initiate the solution
[T]=initiate(n,m,t0);

%Set boundry condition for first time
[T]=Bcs(n,m,T,T0,T1);
Told=T;
error=zeros(MIT);
hold on
%Begin Iteration
while((IT<=MIT)&&(err>eps))
    IT=IT+1
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
    errT=max(max(abs((T-Told))))
    error(IT)=errT;
end

%Postprocces
hold on
subplot(2,1,1)
plot(1:MIT,log10(error),'- r')
xlabel('Iteration')
ylabel('Log(error)')

subplot(2,1,2)
[C1,h1] = contourf(T,30);
text_handle = clabel(C1,h1,'manual');
colorbar
xlabel('x')
ylabel('y')

