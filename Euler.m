function [ q ] = Euler( input, Dt, RHS )
%q returns P,U and V at new time
%   input is value of q in current time
q=input+RHS*Dt;
end

