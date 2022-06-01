function  [Chi, Parm, Px] = quant2(Parm,M,~,~,~,Py,Fy)
% function [Chi, Parm, Px] = quant2(Parm,M,Vars,nn,Inc,Py,Fy)
% 4 dim Hilbert space model
% Parm contains parameters
% Vars contains var no.'s
% nn size of tables
% Py contains relative freq data
% Fy contains raw freq data
% Chi is Chi square diff quant vs saturated
% Px contains predicted probabilities

% model assumes
% each pair forms conjunctions of two events
% each conj is incompatible with another conj
%  each pair uses a different basis
% use unitary matrix to choose basis for a pair

% Parms = AI AS AH IS IH SH A I S H

P = Parm;
P = [ P(1) P(6+[1 2]) ;  ...
      P(2) P(6+[1 3]) ;  ...
      P(3) P(6+[1 4]) ;  ...
      P(4) P(6+[2 3]) ;  ...
      P(5) P(6+[2 4]) ;  ...
      P(6) P(6+[3 4]) ]; 

% P = reshape(P,2,6)';
  

dim = size(M{1},2);

Psy = sqrt(1/dim)*ones(dim,1);   % magnitudes

Px = [];

% ai as ah is ih sh

Mu = cell(1,2);

for j=1:6
    Mu{1} = BuildInComp(P(j,:),M{1},1); 
    Mu{2} = BuildInComp(P(j,:),M{2},2);
    T = real(TwoWayQ(Mu{1},Mu{2},Psy,2,2)); 
        % YY YN
        % NY NN
    n = size(T,1).*size(T,2);
    px = reshape(T,n,1);
    Px = cat(1,Px,px) ;
end

eps = 10^-5;  %-5
Px = eps + (1-2*eps)*Px;
Py = eps + (1-2*eps)*Py;
Chi2 = Fy'*log(Px) ;
Chi1 = Fy'*log(Py) ;
Chi = 2*(Chi1 - Chi2);


