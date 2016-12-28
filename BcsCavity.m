function [P,U,V] = BcsCavity( n,m,P,U,V )
%This function sets Bounry conditions

%Left
U(:,1)=0; %Inlet velocity x component
V(:,1)=0; %Inlet velocity y component
P(:,1)=P(:,2);  %Pressure is interpolated from inside

%Right

U(:,n)=0; %Inlet velocity x component
V(:,n)=0; %Inlet velocity y component
P(:,n)=P(:,n-1);  %Pressure is interpolated from inside

%   Top Wall u=v=0 P(1,j)=P(2,j)
U(1,:)=0; %Inlet velocity x component
V(1,:)=0; %Inlet velocity y component
P(1,:)=P(2,:);  %Pressure is interpolated from inside

%   Bottom Wall u=v=0 P(n,j)=P(n-1,j)
U(m,:)=1; %Inlet velocity x component
V(m,:)=0; %Inlet velocity y component
P(m,:)=P(m-1,:);  %Pressure is interpolated from inside
end

