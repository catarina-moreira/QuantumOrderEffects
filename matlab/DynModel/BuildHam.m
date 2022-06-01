function [ H ] = BuildHam( P,m )
% function [ H ] = BuildHam( P,m )

% m = Hdim
A1 = P(1);
A2 = P(2);
% t=1;

% Hamiltonians
a1 = A1*(1:m)'/m;             % 
a2 = A2*ones((m-1),1) ;

H = diag(a1);
A = diag(a2);
H(1:(m-1),2:m) = H(1:(m-1),2:m) + A;     
H(2:m,1:(m-1)) = H(2:m,1:(m-1)) + A;    % in col out row

end

