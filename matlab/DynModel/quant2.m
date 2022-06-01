function [Chi, parm, Px] = quant2(parm,M,Vars,nn,~,Py,Fy)
% function [Chi, parm, Px] = quant2(parm,M,~,nn,~,Py,Fy)
% change comment lines 73,75 when model fitting to save time

Px = [];   % Predicted
nt = size(Vars,2);

% parameters
A1 = parm(1);   % question A potential
A2 = parm(2);     % diffusion A
% wd = parm(3);  % width of initial state
 wd = 3;   % wd = 3
% t1 = parm(4);  % .5 time int 
% t2 = parm(5);  % 1.0 time int
% t3 = parm(6);  % 1.5 time int
% t4 = parm(7);  % 2.0 time int

% a = parm(4);
a = 1;
t1 = a*(.5); 
t2 = a*(1.0);
t3 = a*(1.5);
t4 = a*(2.0);
% 


tt(1,:) = [t1 t2];  % time intervals .5 and 1
tt(2,:) = [t1 t4];  % time intervals .5 and 2
tt(3,:) = [t3 t2];  % time intervals 1.5 and 1

%    pot diff time
P = [A1   A2   0];

% response scale model
m = size(M{1},2);       % total number of states, make this an odd number
mid = (m+1)/2;          % middle of states

% initial state
% x0 = zeros(m,1);
% x0((mid-wd):(mid+wd),1) = 1;    % initial state
x0 = ((1:m) - mid)';
x0 = (x0./wd).^2;       % normal dist around mid with std = wd
x0 = exp(-x0);
Psy0 = x0./sqrt(x0'*x0);

H = BuildHam(P,m);

for j=1:nt
    v = nn{j};
    ti = tt(j,:);
    U1 = expm(-1i*ti(1)*H) ;  
    Psy1 = U1*Psy0;
    P(3) = -ti(2); 
    % change sign because BuildInComp uses U*M2*U' here we need U'*M2*U
    M2 = BuildInComp(P,M{2});
    T = real(TwoWayQ(M{1},M2,Psy1,v(1),v(2))); % 10 x 10
    n = size(T,1).*size(T,2);
    px = reshape(T,n,1);
    Px = cat(1,Px,px) ;
end


eps = 10^-5;   % this can change results slightly if pred near zero
Px = eps + (1-2*eps)*Px;
Py = eps + (1-2*eps)*Py;

Chi2 = Fy'*log(Px) ;
Chi1 = Fy'*log(Py) ;
Chi = 2*(Chi1 - Chi2);

