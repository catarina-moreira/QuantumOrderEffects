% Main HSM program

clear
clc

% choose data to fit
% data=1 for quant GenerateData , data=2 for TestJoint,
% data = 3 for PSA harm, data = 4 for PSA death

data = 4;

% choose model to fit
model = 2;   %  1= Joint prob, 2 =quant

% choose to fit data or to check fit
fit =3;  % 1 = use particle, 2 = use fminun  3 = check fit

% number of replications of starting points (especially needed for fminunc)
if fit == 3
    nr = 1;
else
    nr = 10;
end



% This example has 4 Attributes A, H , I, U
% some one- way , all pairs, and a couple of reverse orders
% { A, H, I, U, AH, AI, AU, HI, HU, HA, UI}


% your input needed
% ************************
%  enter your # of values for each variable
na=2;  % input # values of A
nh=2;  % input # values of H
ni=2;  % input # values of I
nu=2;  % input # values of U
nv = [na nh ni nu];
Jdim = prod(nv);

% your input needed
% ******************
% read in your raw frequency tables into a cell called Tbls

if data == 1
    % simulated data from quantum model
    load HSdata  %  generated by run GenerateTables
    % loads Tbls
    % 2 x 2 tables have cells,
    % YY YN
    % NY NN
    
    % Data = 1 generates eight 2 x 2 tables, and four one-way tables
    % It generates a cell called Tbls
    % Tbls = { A, H, I, U, AH, AI, AU, HI, HU, HA, UI};
    % A is a 1 x 2 table, B is a 1 x 2 table , ect.
    % A = [fy fn];   fy = freq yes
    % AH is a 2 x 2 table, AI is a 2 x 2 table, ect.
    % AH = [ fyy fyn
    %        fny fnn ] ;
    % the row index is for the variable A, which came first
    % the column index is for variable H, which came second
    
    % your input needed
    % ********************
    % Define variables in tables using Vars
    
    % assign numbers to variables
    % A = 1, H = 2, I = 3, U = 4, 5 = null (for one-way table)
    % The Vars information must be entered to describe variables in tables
    %         null A   null H  null I   null U  A H     A I     A U
    % Vars = { [5 1 ], [5 2 ], [5 3 ], [5 4 ], [1 2 ], [1 3 ], [1 4 ], ...
    %          [2 3 ], [2 4 ], [3 4 ], [2 1 ], [4 3 ] };
    %           H I     H U     I U     H A     U I
    % Vars is included in HSdata
    
    % define compatibility relations
    
    % var 1 is comp with var 3 , forms initial basis
    % var 1 is first cubit, var 3 is second
    % var 1 is rotated to var 2
    % var 3 is rotated to var 4
    %     1st 2nd cubit
    % var 1 rotated to 2, var 3 rotated to 4
    Comp = [1 3];
    Inc =  [2 4];
    
    
    % This example uses arbitrary complex initial state
    % 3 magnitudes, 3 phases  for a 4-dim initial state
    %  one mag arbitrary becasue length = 1
    %  one phase arbitrary because prob's do not depend on common phase
    % Arbitrary complex Hemitian matrices in off diag
    % each 2 x 2 Ham has 3 parms, 2 Ham x 3 parms = 6
    % total no. of parameters  [ six initial, six Ham]
    np = 12;
    
    
elseif data == 2
    % simulated data from joint probability model
    % load saved data tables, called NTbls, using a cell for each table
    % loads data for which the joint prob is true, to test program
    load TestJoint
    Tbls = N2Tbls;
    % same type of tables as above but generated from Joint
    % contains NTbls, Px, Vars, nn
    % Px is the joint probability distribution
    % NTbls is the cell with Tables
    
elseif data == 3
    % real data ratings of PSA's by 184 people
    % Harm PSA's
    load PSAharmO1
    % loads Tbls cell with 6 cells
    % 4 Attributes:  Persuasive, Informative, Believable, Likeable,
    % 6, 2 by 2 tables
    % Tbls = {pi pb pl ib il bl }
    % 1 = p, 2 = i, 3 = b, 4 =l
    %           1     2     3     4     5     6
    % Vars = {[1 2],[1 3],[1 4],[2 3],[2 4],[3 4]}
    % Vars is included in PSAH
    % 8 obs per person, 184 participants, N=1472 observations per table
    
    
    % define compatibility relations  1=p 2=i 3=b 4=l
    % var 2 is comp with var 3 , forms initial basis
    % var 3 is rotated to var 1
    % var 2 is rotated to var 4
    Comp = [2 3]; % [2 3];
    Inc =  [4 1]; % [4 1];
    
    % This example uses real valued initial state
    % The initial state is set equal to one of the 2 by 2 tables
    % The Hamiltonian is real
    % one 2 x 2 Ham for each inc var
    % 2 real parm's used for each Ham
    %  2 x 2 = 4 Ham parms
    % no. of parameters
    np = 4;
    
else
    % real data ratings of PSA's by 184 people
    % Harm PSA's
    load PSADeathO1
    % loads Tbls cell with 6 cells
    % 4 Attributes:  Persuasive, Informative, Believable, Likeable,
    % 6, 2 by 2 tables
    % Tbls = {pi pb pl ib il bl }
    % 1 = p, 2 = i, 3 = b, 4 =l
    %           1     2     3     4     5     6
    % Vars = {[1 2],[1 3],[1 4],[2 3],[2 4],[3 4]}
    % Vars is included in PSAH
    % 8 obs per person, 184 participants, N=1472 observations per table
    
    
    % define compatibility relations
    % var 2 is comp with var 3 , forms initial basis
    % var 3 is rotated to var 1
    % var 2 is rotated to var 4
    Comp = [2 3]; % [2 3];
    Inc =  [4 1]; % [4 1];
    
    % This example uses real valued initial state
    % The initial state is set equal to one of the 2 by 2 tables
    % The Hamiltonian is real
    % one 2 x 2 Ham for each inc var
    % 2 real parm's used for each Ham
    %  2 x 2 = 4 Ham parms
    % no. of parameters
    np = 4;   % 4 real
    
end


% **********************
% Define size of Tbls using nn

nt = size(Tbls,2);   % no. of tables
nn = cell(1,nt);
for j = 1:nt
    nn{j} = size(Tbls{j});  % table size
end


% Next step creates a vector of data to be predicted

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

if  model == 1   % test joint probability model
    
    M4 = ProjJP(nv);
    
    lb = zeros(Jdim,1);
    ub = ones(Jdim,1);
    options = optimoptions('fmincon','MaxFunEvals',10^10,'MaxIter',1000,'Display','off','Algorithm','sqp');
    parm0 = rand(Jdim,1); parm0 = parm0./sum(parm0);
    [Parm, Chi] = fmincon(@(parm) JointP(parm,M4,Vars,nn,Py,Fy),parm0,[],[],ones(1,Jdim),1,lb,ub,[],options);
    
    df = (ny -nt) - (Jdim-1);
    p = 1-cdf('chi2',Chi,df);
    disp('Joint Prob Model')
    disp('Chisq diff       df       p')
    disp([Chi df p])
    
elseif model == 2    % quantum model with 12 parameters used in articl
    
    M2 = ProjQP(nv,Comp,Inc,[]);
    
    lb = -5*ones(np,1) ;  % only for swarm
    ub =  5*ones(np,1) ;
    
    %  use for particle swarm
    fun = @(parm) quant2(parm,M2,Vars,nn,Comp,Inc,Py,Fy);
    options1 = optimoptions('particleswarm','SwarmSize',100,'UseParallel',true,'Display','iter','MaxIter',200);
    
    % use for fminunc
    options2 = optimoptions('fminunc','MaxFunEvals',10^10,'MaxIter',1000,'Display','off','Algorithm','quasi-newton');
    
    
    PR = zeros(nr,np); ChiR = zeros(nr,1);
    
    for rep = 1:nr
        
        if fit ==1
            [PR(rep,:), ChiR(rep), Px] = particleswarm(fun,np,lb,ub,options1);
            
        elseif fit == 2
            % HSdata
            %                parm0 = [ -0.2168   -0.5833    0.2752   2.2920   0.9383   0.0400 ...
            %                      -0.5911   -0.5037    0.8862 ...
            %                      -1.2405  -0.4334     1.2976 ];   % best fit parms HSD data
            
            parm0 = zeros(1,np);
            parm0 = parm0 + 1*randn(1,np);
            
            [PR(rep,:), ChiR(rep)] = fminunc(@(parm) quant2(parm,M2,Vars,nn,Comp,Inc,Py,Fy),parm0,options2);
            [rep ChiR(rep)]
            
        else
            if data == 1
                parm = [  -0.2168   -0.5833    0.2752   2.2920   0.9383   0.0400 ...
                    -0.5911   -0.5037    0.8862 ...
                    -1.2405  -0.4334     1.2976];   % best fit parms HSData Chi = 9.3346e-05
            elseif data ==3  % harm
                % PSA data harm Chi = 171 using all 12 below
                %                 parm = [ -0.8996 -0.7351 -0.0101 0.8458 0.1841 0.7529 2.5094 -2.1927 -0.2803 1.1012 -0.3855 0.1527];
                parm = [ -0.2106    0.5388   -0.4769   -0.3580];  % chi= 194 using 4
                
            else  % death
                % chi = 145.66 using all 12 parms
                %  parm = [-0.9973 -0.7592 -0.3713 -0.8411 -3.8289 1.7095 -4.9620 3.0067 3.2081 -2.1861 0.4861 -3.1524 ] ;
                parm = [-0.0614   -0.5502    0.6993    0.4028]; %  Chi = 150
                
            end
            PR = parm;
            
        end  % pick fit method
        
    end % rep
    
    [Chi, Ind] = min(ChiR);    % pick best fit
    parm = PR(Ind,:);
    
    [Chi,parm, Px] = quant2(parm,M2,Vars,nn,Comp,Inc,Py,Fy);
    
    disp('Quantum model')
    disp('Chisquare difference')
    disp(Chi)
    if data == 1
        disp('parameters for initial state')
        disp(parm(1:6))
        disp('parameters for each inc var Ham')
        disp(parm(7:12))
    else
        disp('parameters for each inc var Ham')
        disp(parm)    
    end
    
    'predicted left ; observed right'
    % each table stretched out rows for 1st col then rows for second col 
    %  in same order of Tbls read into PY 
    % e.g. for first 2 x 2 table in Data == 4
    %  [781 256
    %    94 341]
    % stretched out to
    %  [ 781
    %     94
    %    256
    %    341 ]
    % see line 195 
    [Px Py]
    
end % model



