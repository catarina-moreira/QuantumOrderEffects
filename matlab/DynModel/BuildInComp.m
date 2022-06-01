function [ M ] = BuildInComp( P,X)
% function [ M ] = BuildInComp( P,X )



% X contains M for inc;

%  m = size(c,2);   % no. values of C
%  n = max(c{m});   % dim of H

m = size(X,1);
n = size(X,2);

H = BuildHam(P,n);
U = expm(-1i*P(3)*H);

C = zeros(m,n,n) ;
for j = 1:m
    % C(j,:,:) = X(j,:,:);
    C(j,:,:) = U*squeeze(X(j,:,:))*U';
end

M = C;
end

