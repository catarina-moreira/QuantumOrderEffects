function [ M ] = BuildInComp( P,X )
% function [ M ] = BuildInComp( P,X )

%   |A> * |B>  compatible space
%   na = dim of A
%   nb = dim of B
%   c = subspaces for measurements of c, row of cells

% X contains {M1,M2,M3,1};

%  m = size(c,2);   % no. values of C
%  n = max(c{m});   % dim of H

% t =  P(5);
t= 1;

m = size(X,1);
n = size(X,2);

[Ha, Hb] = BuildHam(P);


Ua = expm(-1i*t*Ha);
Ub = expm(-1i*t*Hb);

U = kron(Ua,Ub);



C = zeros(m,n,n) ;
for j = 1:m
    % C(j,:,:) = X(j,:,:);
    C(j,:,:) = U*squeeze(X(j,:,:))*U';
end

M = C;
end

