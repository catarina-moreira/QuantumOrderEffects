function K = BuildInt(P,m)
% K = BuildInt(P,m)

 a = P(1); % off diag top
 b = P(2); % off diag bottom 

 
gam = -(a+b)*ones(m,1);
alfa = a*ones(m,1);    % rate up = (prob up)/h
beta = b*ones(m,1);    % rate dn = (prob dn)/h

up = alfa(2:m,1)  ;
dn = beta(1:(m-1),1);
cn = gam.*ones(m,1);

A = diag(up);
B = diag(dn);
K = diag(cn);
K(1:(m-1),2:m) = K(1:(m-1),2:m) + A;     
K(2:m,1:(m-1)) = K(2:m,1:(m-1)) + B;    
K(1,1) = -K(2,1);
K(m,m) = -K(m-1,m);

% in col, out row
% col's sum to zero
