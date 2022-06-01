function  [Chi, Parm, Px] = quant2(Parm,M,Vars,nn,Inc,Py,Fy)
% function [Chi, Parm, Px] = quant2(Parm,M,Vars,nn,Inc,Py,Fy)
% six dim Hilbert space model
% P contains parameters
%  first Dim-1 parms for initial, one coordinate set = 1
% M contains projectors
% Vars list of var's
% nn contains sizes of tables
% Inc indicates incompatible vars
% Py  contains relative freq vector
% Fy contains raw freq vector
% Chi is Chi square diff quant vs saturated
% Px contains predicted probabilities


% Assign parameters
Px = [];

np = size(Parm,2);
nt = size(nn,2);
ncat = nn{nt}(1);
ndec = nn{nt}(2);

P1 = [ 1 Parm(1:(ncat-1)) ]' ;

Psy = kron(P1,ones(ndec,1));
Psy = Psy./sqrt(Psy'*Psy);

nt = size(Vars,2);

M{Inc} = BuildInComp(Parm(2:np),M{Inc});

for j=1:nt 
    rc = Vars{j};
    v = nn{j};
    T = real(TwoWayQ(M{rc(1)},M{rc(2)},Psy,v(1),v(2))) ;  
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


