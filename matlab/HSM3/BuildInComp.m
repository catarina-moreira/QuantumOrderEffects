function [ M ] = BuildInComp( P,X, pos )
% function [ M ] = BuildInComp( P,X,pos )

%   |A> * |B>  compatible space
%   na = dim of A
%   nb = dim of B

m = size(X,1);
n = size(X,2);

H = BuildHam(P,pos); % pos = 1st or 2nd
U = expm(-1i*(pi/2)*H);

C = zeros(m,n,n) ;
for j = 1:m
    % C(j,:,:) = X(j,:,:);
    C(j,:,:) = U*squeeze(X(j,:,:))*U';
end

M = C;

end

