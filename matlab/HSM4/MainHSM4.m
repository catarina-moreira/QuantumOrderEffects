% Main HSM program
% Application to A , B,  AB tables
% 2 levels A
% 2 levels B
% Example using Face data from Busemeyer, Wang, Mogiliansky, JMP 2009

clear
clc


% choose model to fit
model = 2;   %  1= Joint prob, 2 =quant

% choose to fit data or to check fit
fit = 2 ;  % 0 = compute predictions only, 1 = use particle, 2 = fminun
% number of replications of starting points (especially needed for fminunc)
nr = (fit>0).*5 + (fit==0).*1;

% This example has 2 Attributes A, B

% your input needed
% ******************************
%  enter your # of values for each varFiable
%  set extra vars equal to one
na=2;  % input # values of A
nb=2;  % input # values of B
nc=1;  % input # values of C  not used
nd=1;  % input # values of D  not used
nv = [na nb nc nd];
Jdim = prod(nv);


% your input needed
% ******************
% read in your tables into a cell called Tbls

load FaceData

ch = 1;  % ch = 1 bad, ch = 2 good

% TblsG  TblsG = {TCg, TDg, TCDg} ;    % fit good faces
% TBlsB  TblsB = {TCb, TDb, TCDb} ;    % fit bad faces

if ch == 1
    Tbls = TblsB;
else
    Tbls = TblsG;
end

% cat - dec results from JMP 2009
% TDb  F D   decision alone bad faces
% TDg  F D   decision alone good faces
% TCb  G B   cat alone bad faces
% TCg  G B   cat alone good faces

% TCDb      G   B   % C&D  bad faces
%  GF GD BF  BD
%  4   1  3  10

%         F   D
%    G   4   1
%    B   3   10



% your input needed
% ********************
% Define variables in tables using Vars
% numeric coding for variable names


% 1 = A for tables including A   Categ
% 2 = B for tables including B   Action
% 3 = C for tables including C   Change Basis for Action
% 4 = null for one way table
%  e.g., CA is 3 for C, 1 for A
%         null Dec  Cat Dec
% Vars = { [4 3], [1 3]};
Vars = { [4 1], [4 3], [1 3]};

% define size of tables using nn
nt = size(Tbls,2);   % no. of tables
nn = cell(1,nt);
for j = 1:nt
    nn{j} = size(Tbls{j});
end

%  this example used proportions, they needed to be changed to freq
Ns =  1;
% Ns = 1 if you enter frequencies to begin
% Ns = sample size if you enter proportions
% Ns>1 only if data are proportions rather than freq

% create vector of data to be predicted

% Fy
% [ G B F D FG FB DG DB]

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

if  model == 1   % test joint probability model for a single type face
    
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
    
elseif model == 2    % quantum model with 7 parameters
    Comp = [1 2];
    Inc = 3; % define incompatible relations
    M2 = ProjQP(nv,Comp,Inc,[]);
    
    np = 5;
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
            
            parm0 = [-1.8545    1.9717    4.0830   -4.9978   -3.6899];
            parm0 = parm0 + 1*randn(1,np);
            [PR(rep,:), ChiR(rep)] = fminunc(@(parm) quant2(parm,M2,Vars,nn,Inc,Py,Fy),parm0,options2);
            [rep ChiR(rep)]
            
        else
            if ch ==1
                PR(rep,:) = [  -1.8545    2.8874   -4.9315    1.3223    1.0980];  % bad face
                ChiR(rep) =  0.0515 ;
            else
                PR(rep,:) = [ -0.5207   -3.3433    4.9474   -1.2647    1.4689 ];  % good face
                ChiR(rep) =  0.0450 ;
            end
            
            
            
        end  % pick fit method
        
    end % rep
    
    [Chi, Ind] = min(ChiR);    % pick best fit
    parm = PR(Ind,:);
    
    [Chi,parm, Px] = quant2(parm,M2,Vars,nn,Inc,Py,Fy);
    
    
    disp('Quantum model')
    disp('Chisquare difference')
    disp(Chi)
    disp('parameters')
    disp(parm)
    'predicted left observed right'
    ' [ F D FG FB DG DB] '
    
    [Px Py]
    
end % model



