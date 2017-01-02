function [ disspx ] = DisspX( i,j,m,Q,epsx )
%this function calculates diisspation Q in X direction for each point for
%incommpresible flow
% Set an intermediate component to zero if cell is cell adjucent to the
% boundry
if (j==2)  
    qijm2=0;  %set Q(i,j-2)=0
else
    qijm2=Q(i,j-2);
end
if (j==m-1)  
    qijp2=0;  %set Q(i,j+2)=0
else
    qijp2=Q(i,j+2); 
end
disspx=epsx*( qijm2-4*Q(i,j+1) + 6*Q(i,j)-4*Q(i,j-1)+ qijp2 );

end

