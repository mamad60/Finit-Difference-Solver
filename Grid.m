function [ X,Y,dL,dH ] = Grid(m,n,L,H )
%This  function compute location of each grid points
%n,m no. of nodes in x,y direction
%Lenght and Height of the Channel


dL=L/(m-1);
dH=H/(n-1);

for j=1:m
    for i=1:n
        X(i,j)=(j-1)*dL;
        Y(i,j)=(i-1)*dH;
    end
end
end   


