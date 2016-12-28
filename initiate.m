function [P,U,V] = initiate(n,m,p0,u0,v0)
for i=1:n
    for j=1:m
       P(i,j)=p0;
       U(i,j)=u0;
       V(i,j)=v0;
    end
end
end

