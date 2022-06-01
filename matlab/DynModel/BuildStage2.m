function [ M ] = BuildStage2( P,X)
% function [ M ] = BuildStage2( P,X )



% X contains M for inc;

%  m = size(c,2);   % no. values of C
%  n = max(c{m});   % dim of H

m = size(X,1);
n = size(X,2);

K = BuildInt(P,n);
T = expm(P(3)*K);

C = zeros(m,n,n) ;
for j = 1:m
    % C(j,:,:) = X(j,:,:);
    C(j,:,:) = squeeze(X(j,:,:))*T;
end

M = C;
end

