function [ M ] = BuildInComp( P,X )
% function [ M ] = BuildInComp( P,X )

%   |A> * |B>  compatible space
%   na = dim of A
%   nb = dim of B
%   c = subspaces for measurements of c, row of cells

% X contains {M1,M2,M3,1};

%  m = size(c,2);   % no. values of C
%  n = max(c{m});   % dim of H

t =  pi/2;

m = size(X,1);
n = size(X,2);

H = BuildHam(P);
U = expm(-1i*t*H);

C = zeros(m,n,n) ;
for j = 1:m
    % C(j,:,:) = X(j,:,:);
    C(j,:,:) = U*squeeze(X(j,:,:))*U';
end

M = C;
end

