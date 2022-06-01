% Main HSM program

clear
clc

% choose data to fit
% select data set  data = 1 CSdata, data = 2 joint prob

data = 1;

% choose model to fit
model = 2;   %  1= Joint prob, 2 =quant

% choose to fit data or to check fit
fit = 3;  % 1 = use particle, 2 = use fminun  3 = check fit
% number of replications of starting points (especially needed for fminunc)
nr = 1;

% This example has 3 Attributes A, B , C

% your input needed 
% ******************************
%  enter your # of values for each variable
%  set extra vars equal to one 
na=2;  % input # values of A
nb=3;  % input # values of B
nc=4;  % input # values of C
nd=1;  % input # values of D  
nv = [na nb nc nd];
Jdim = prod(nv);


% your input needed
% ******************
% read in your tables into a cell called Tbls

if data == 1   % data generate by     run CSDesign
    load CSdata
    
    % Tbls containing 6 Tables of frequencies
    % TC =  real(OneWay(MC,Psy,nc))           1 x 4
    % TAB = real(TwoWay(MA,MB,Psy,na,nb))     2 x 3
    % TAC = real(TwoWay(MA,MC,Psy,na,nc))     2 x 4
    % TCA = real(TwoWay(MC,MA,Psy,nc,na))     4 x 2
    % TBC = real(TwoWay(MB,MC,Psy,nb,nc))     3 x 4
    % TCB = real(TwoWay(MC,MB,Psy,nc,nb))     4 x 3
    % Tbls = {TC,  TAB,      TAC,   TCA',    TBC,     TCB };
    
else 
    % load saved data tables, called NTbls, using a cell for each table
    load TestJoint
    Tbls = N2Tbls;
    % same type of tables as CSDesign but generated from Joint
    % contains NTbls, Px, Vars, nn

end

% your input needed
% ********************
% Define variables in tables using Vars
% numeric coding for variable names

% 1 = A for tables including A
% 2 = B for tables including B
% 3 = C for tables including C
% 4 = null for one way table
%  third index indicates first vs second order
%  e.g., table C A is 3 for C, 1 for A,  producing [3 1]
%  e.g., table A C is 1 for A, 3 for C,  producing [1 3]
%          null C    A B      A C      C A       B C     C B
% Vars = { [4  3], [1 2],    [1 3 ],  [3 1],   [2 3],   [3 2]};


% define size of tables using nn 
nt = size(Tbls,2);   % no. of tables
nn = cell(1,nt);
for j = 1:nt
    nn{j} = size(Tbls{j});
end

%  this example used proportions, they needed to be changed to freq
Ns =  100; 
% Ns = 1 if you enter frequencies to begin 
% Ns = sample size if you enter proportions
% Ns>1 only if data are proportions rather than freq

% create vector of data to be predicted
Py = []; Fy = [];
for j = 1:nt
    T = Tbls{j};
    n = size(T,1).*size(T,2);
    fy = Ns*reshape(T,n,1);
    py = fy./sum(fy);
    Fy = cat(1,Fy,fy);
    Py = cat(1,Py,py) ;
end
ny = size(Py,1);

if  model == 1   % test joint probability model
    
    M3 = ProjJP(nv);
    
    lb = zeros(Jdim,1);
    ub = ones(Jdim,1);
    options = optimoptions('fmincon','MaxFunEvals',10^10,'MaxIter',1000,'Display','off','Algorithm','sqp');
    parm0 = rand(Jdim,1); parm0 = parm0./sum(parm0);
    [Parm, Chi] = fmincon(@(parm) JointP(parm,M3,Vars,nn,Py,Fy),parm0,[],[],ones(1,Jdim),1,lb,ub,[],options);
    
    df = (ny -nt) - (Jdim-1);
    p = 1-cdf('chi2',Chi,df);
    disp('Joint Prob Model')
    disp('Chisq diff       df       p')
    disp([Chi df p])
    
elseif model == 2    % quantum model with 12 parameters used in article
    Comp = [1 2];
    Inc = 3; % define incompatible relations
    M2 = ProjQP(nv,Comp,Inc,[]);
    
    np = 11;
    lb = -5*ones(np,1) ;
    ub =  5*ones(np,1) ;
    
    %  use for particle swarm
    fun = @(parm) quant2(parm,M2,Vars,nn,Inc,Py,Fy);
    options1 = optimoptions('particleswarm','SwarmSize',100,'UseParallel',true,'Display','iter','MaxIter',200);
    
    % use for fminunc
    options2 = optimoptions('fminunc','MaxFunEvals',10^10,'MaxIter',1000,'Display','off','Algorithm','quasi-newton');
   
    PR = zeros(nr,np); ChiR = zeros(nr,1);
    
    for rep = 1:nr
        
        if fit ==1
            [PR(rep,:), ChiR(rep)] = particleswarm(fun,np,lb,ub,options1);
            
        elseif fit == 2
            P1 = [    0.2685 -0.7600 0.0883 0.3514 -0.1935 -0.4260]; 
            ma = 1/3; sa = 1; mb = -1/15; sb = 2.1; t = pi/2;
            P2 = [ma sa mb sb t];
            parm0 = [P1 P2];
            parm0 = parm0 + .2*randn(1,11);
            [PR(rep,:), ChiR(rep)] = fminunc(@(parm) quant2(parm,M2,Vars,nn,Inc,Py,Fy),parm0,options2);
            [rep ChiR(rep)]
            
        else
            if data == 1
            P1 = [    0.2685 -0.7600 0.0883 0.3514 -0.1935 -0.4260]; 
          %  ma = 1/3; sa = 1; mb = -1/15; sb = 2.1; 
            t = pi/2;
            ma = t*(1/3); sa = t*1; mb = t*(-1/15); sb = t*(2.1);

            P2 = [ma sa mb sb];
            PR = [P1 P2];
            end
            
            
            
        end  % pick fit method
        
    end % rep
    
    [Chi, Ind] = min(ChiR);    % pick best fit
    parm = PR(Ind,:);
    
                [Chi,parm, Px] = quant2(parm,M2,Vars,nn,Inc,Py,Fy);

    
    disp('Quantum model')
    disp('Chisquare difference')
    disp(Chi)
    disp('parameters for initial state')
    % note that in the program they are normalized
    % so the normalized initial state equals
    x = parm(1:6)'; x = x./sqrt(x'*x);
    disp(x)
    disp('parameters for Hamiltonians')
    P = parm(7:10);
    disp(P)
    % first 2 are for HA
    % second 2 are for HB
    
disp('Transitions from A to C')
abs(Ua).^2
disp('Transitions from B to C')
abs(Ub).^2
    
      'predicted left observed right'
    % each table stretched out rows for 1st col then rows for second col 
    %  in same order of Tbls read into PY 
    % e.g. for second 2 x 3 table in Data == 1 
   %  0.0721    0.5777    0.0078
   %  0.1235    0.0374    0.1815
    % stretched out to
   % 0.0721
   % 0.1235
   % 0.5777
   % 0.0374
   % 0.0078
   % 0.1815
   % see line 95
      [Px Py]
    
end % model



