% Quantum question order model for rating scalse
%
% In this example
% Fit SO OS distributions for 9 point rating scales
% 99 states
% uses inverse to transform back to neutral
% Random starting point
% Five parameters
clear
clc

%  enter your # of values for each variable
ns=9;  % input # values of self ratings
no=9;  % input # values of other ratings
nc=1;  % input #   fix to one 
nd=1;  % input #   fix to one 
nv = [ns no nc nd];
Jdim = prod(nv);
Hdim = 99;

% JointCoact contains two joint freq tables, 
% one for AB order, other for BA order
% this example defines FA for the AB table, FB for the BA table
 
load JointCoacT   
% this file contains a cell Tbls with two tables
% FS % self  Tbls{1}
% FO % other Tbls{2}


% FA self rows, other columns, self first 
% FB other rows, self columns, other first
% 1 = self,  2 = other
% Vars = { [1 2], [2 1] }

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define size of tables using nn 
nt = size(Tbls,2);   % no. of tables
nn = cell(1,nt);
for j = 1:nt
    nn{j} = size(Tbls{j});
end


% create vector of data to be predicted
Py = []; Fy = [];
for j = 1:nt
    T = Tbls{j};
    n = size(T,1).*size(T,2);
    fy = reshape(T,n,1);
    py = fy./sum(fy);
    Fy = cat(1,Fy,fy);
    Py = cat(1,Py,py) ;
end
ny = size(Py,1);

Comp = 0;  % both incompatible with initial state
Inc = [1 2];
M2 = ProjQP(nv,Comp,Inc,Hdim);


% options = optimset('Display','iter','MaxFunEvals',2000,'MaxIter',1000,'TolX',1e-6);
  options = optimset('Display','off','MaxFunEvals',2000,'MaxIter',1000,'TolX',1e-3);


% initial start values
%  b1    b2      a1         a2    t=1

% Coactive

 % min SSE
 
 
 
 % N = 133
 % SSE
 
 parm0 = [ -14.5725  -16.7430    4.5975    4.4946    2.7558 ]; % R2 = .9046, SSE= .0086
 % Chi 
% parm0 = [-9.5695  -13.3291    4.6880    4.4575    2.2297 ] ;  % Chi = 854.4356
 
% parm0 = [ -13.2794  -17.5574    4.7413    4.5264    2.2768];  % Chi = 839.2622

 
sw = 0;
% sw = 0 turn off fitting routine,
% sw = 1 turn on fitting routing

ChiV = [];
ParmM = [];

reps = 1;   % number of random starting positions
std = sw*(0);    % amount of jitter

for n = 1:reps
    nj = size(parm0,2);
    jitter = std*randn(1,nj);
        
    parm0 = parm0 + jitter;
    
    if sw == 1
        
        [parm,Chi] = fminsearch(@(parm) ABmodel(parm,M2,Vars,nn,Inc,Py,Fy), parm0, options);
        [ n Chi parm]
        ChiV = [ChiV ; Chi];
        ParmM = [ParmM ; parm];
    else
        parm = parm0;
    end
end

[Chi,Px] = ABmodel(parm,M2,Vars,nn,Inc,Py,Fy);
Chi

% save QFitResults parm PM Fdata

disp('Parms')
A1 = parm(1);   % question A potential
B1 = parm(2);   % question B potential
A2 = exp(parm(3));     % diffusion A
B2 = exp(parm(4));     % diffusion B
mix = 1/(1+exp(-parm(5)));
disp('   A-potent  B-potent  A-drift   B-Drift    Mix')
disp([A1 B1 A2 B2 mix])

'predicted '
Px = reshape(Px,9,9,2)  % matches old predictions (except other not transposed)
'observed'
Py = reshape(Py,9,9,2);
