function [T] = Bcs( n,m,T,T0,T1)
%This function sets Bounry conditions

%Left
T(:,1)=T0; %Dirichlet Bcs On Left side

%Right
T(:,m)=T1; %Dirichlet Bcs On Right side

%Top
T(1,:)=T(2,:); %Neumann Bcs On Top side

%Bottom
T(n,:)=T(n-1,:); %Neumann Bcs Bcs On Bottom side
end

