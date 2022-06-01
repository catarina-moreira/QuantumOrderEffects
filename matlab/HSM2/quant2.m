function  [Chi, Parm, Px] = quant2(Parm,M,Vars,nn,Inc,Py,Fy)
% function [Chi, Parm, Px] = quant2(Parm,M,Vars,nn,Inc,Py,Fy)
% six dim Hilbert space model
% P contains parameters
% M contains projectors
% Vars list of var's
% nn contains sizes of tables
% Inc indicates incompatible vars
% Py  contains relative freq vector
% Fy contains raw freq vector
% Chi is Chi square diff quant vs saturated
% Px contains predicted probabilities

% model assumes
% A=1 and B=2 compatible
% C=3 incompatible
% C=3 Inc A C=3 Inc B;

% Assign parameters
P = Parm;
Px = [];

np = size(P,2);
dim = size(M{1},2);
Psy = P(1,1:dim)'; Psy = Psy./sqrt(Psy'*Psy);
nt = size(Vars,2);
rp = dim + (1:(np-dim));

M{Inc} = BuildInComp(P(1,rp),M{Inc});

for j=1:nt  
    rc = Vars{j};
    v = nn{j};
    T = real(TwoWayQ(M{rc(1)},M{rc(2)},Psy,v(1),v(2)));   % 2 x 3
    n = size(T,1).*size(T,2);
    px = reshape(T,n,1);
    Px = cat(1,Px,px) ;
end

eps = 10^-5;
Px = eps + (1-2*eps)*Px;
Py = eps + (1-2*eps)*Py;
Chi2 = Fy'*log(Px) ;
Chi1 = Fy'*log(Py) ; 
Chi = 2*(Chi1 - Chi2);


