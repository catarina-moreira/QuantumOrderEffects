function [Chi, Px] = ABmodel(parm,M,Vars,nn,Inc,Py,Fy)
% function [Chi, Px] = ABmodel(parm,M,Vars,nn,Inc,Py,Fy)
% change comment lines 73,75 when model fitting to save time

Px = [];   % Predicted
nt = size(Vars,2);

% parameters
A1 = parm(1);   % question A potential
B1 = parm(2);   % question B potential
A2 = exp(parm(3));     % diffusion A
B2 = exp(parm(4));     % diffusion B
mix = 1/(1+exp(-parm(5)));

% response scale model
m = size(M{1},2);       % total number of states, make this an odd number
mid = (m+1)/2;          % middle of states
n2 = size(M{1},1);      % no. rating values

% initial state
wd = 5;
x0 = zeros(m,1);
x0((mid-wd):(mid+wd),1) = 1;    % initial state
Psy = x0./sqrt(x0'*x0);

M{Inc(1)} = BuildInComp([A1 A2],M{Inc(1)});
M{Inc(2)} = BuildInComp([B1 B2],M{Inc(2)});

for j=1:nt
    rc = Vars{j};
    v = nn{j};
    T = real(TwoWayQ(M{rc(1)},M{rc(2)},Psy,v(1),v(2)));   % 9 x 9
    
    % adjustment for bias to respond in the middle of scale
    Mid = (n2+1)/2;
    Mid = [ floor(Mid) ceil(Mid) ];
    M5 = zeros(n2);
    M5(Mid,Mid) = 1/(1+ (Mid(1)-Mid(2)));
    T = mix*T + (1-mix)*M5;
    % end adjustment
    
    n = size(T,1).*size(T,2);
    px = reshape(T,n,1);
    Px = cat(1,Px,px) ;
end

eps = 10^-5;
Px = eps + (1-2*eps)*Px;
Py = eps + (1-2*eps)*Py;

Chi2 = Fy'*log(Px) 
Chi1 = Fy'*log(Py) 
Chi = 2*(Chi1 - Chi2);

