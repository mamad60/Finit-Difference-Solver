function [ disspy ] = DisspY( i,j,n,Q )
%this function calculates diisspation Q in X direction for each point for
%incommpresible flow
% Set an intermediate component to zero if cell is cell adjucent to the
% boundry
if (i==2)  
    qim2j=0;  %set Q(i-2,j)=0
else
    qim2j=Q(i-2,j);
end
if (i==n-1)  
    qip2j=0;  %set Q(i+2,j)=0
else
    qip2j=Q(i+2,j); 
end
epsy=.001; %diisipation coefficient

disspy=epsy*( qim2j- 4*Q(i+1,j) + 6*Q(i,j)-4*Q(i-1,j)+ qip2j );

end

