function [P,U,V] = Bcs( n,m,P,U,V )
%This function sets Boundary conditions

%Inlet
U(:,1)=1; %Inlet velocity x component
V(:,1)=0; %Inlet velocity y component
P(:,1)=2*P(:,2)-P(:,3);  %Pressure is interpolated from inside

%Outlet 
P(:,m)=0;
U(:,m)=2*U(:,m-1)-U(:,m-2);  %U vel is interpolated from inside
V(:,m)=2*V(:,m-1)-V(:,m-2);  %V vel is interpolated from inside

%   Top Wall u=v=0 P(1,j)=P(2,j)
U(1,:)=0;
V(1,:)=0;
P(1,:)=P(2,:);

%   Bottom Wall u=v=0 P(n,j)=P(n-1,j)
U(n,:)=0;
V(n,:)=0;
P(n,:)=P(n-1,:);
end

