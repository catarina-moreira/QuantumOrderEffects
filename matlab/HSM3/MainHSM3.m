% Main HSM program

clear
clc

% choose data to fit
data = 1;

% choose model to fit
model = 2;   %  1= Joint prob, 2 =quant

% choose to fit data or to check fit
fit = 3;  % 1 = use particle, 2 = use fminun  3 = check fit
% number of replications of starting points (especially needed for fminunc)
nr = 1;

% table = 6 is used in this example: female avatar sexualized
% This example has 4 Attributes A, I , S, H
% Attractive Intelligent Social Honest
% 6 pairs of tables, only one order used in this application
% AI AS AH IS IH SH

% your input needed
% *****************************
%  enter your # of values for each variable
na=2;  % input # values of A
ni=2;  % input # values of I
ns=2;  % input # values of S
nh=2;  % input # values of H
nv = [na ni ns nh];
Jdim = prod(nv);


% your input needed
% ******************
% read in your raw frequency tables into a cell called Tbls


if data==1
    
    load AvatarFemSex
    % contains Tbls Vars
    
    % Tbls 6 tables
    % ai as ah is ih sh
    
    % each table is 2 x 2
    % YY YN
    % NY NN
    
elseif data ==2
    
    load AvatarFemNotS
    % contains Tbls Vars
    
    % Tbls 6 tables
    % ai as ah is ih sh
    
    % each table is 2 x 2
    % YY YN
    % NY NN
    
elseif data == 3
    % real data ratings of PSA's by 184 people
    % Harm PSA's
    load PSAdeathO1
    % loads Tbls cell with 6 cells
    % 4 Attributes:  Persuasive, Informative, Believable, Likeable,
    % 6, 2 by 2 tables
    % Tbls = {pi pb pl ib il bl }
    % 1 = p, 2 = i, 3 = b, 4 =l
    %           1     2     3     4     5     6
    % Vars = {[1 2],[1 3],[1 4],[2 3],[2 4],[3 4]}
    % Vars is included in PSAH
    % 8 obs per person, 184 participants, N=1472 observations per table
else
    load PSAharmO1
    % same as PSAdeath
end

% your input needed
% ********************
% Define variables in tables using Vars

% assign numbers to variables
% A = 1, I = 2, S = 3, H = 4, 5 = null (for one-way table)
% The Vars information must be entered to describe variables in tables
%         null A   null I  null S   null H  A I     A S     A H
% Vars = { [5 1 ], [5 2 ], [5 3 ], [5 4 ], [1 2 ], [1 3 ], [1 4 ], ...
%          [2 3 ], [2 4 ], [3 4 ] };
%           I S     I H     S H
% Vars is included in HS3dat



% **********************

% Define size of Tbls using nn

nt = size(Tbls,2);   % no. of tables
nn = cell(1,nt);
for j = 1:nt
    nn{j} = size(Tbls{j});
end


% build data vector to fit
Py = []; Fy = [];
for j = 1:nt
    T = Tbls{j};
    n = size(T,1).*size(T,2);
    fy = reshape(T,n,1);
    py = fy./sum(fy);
    Fy = cat(1,Fy,fy);  % 24 x 1 col: YY NY YN NN YY NY YN NN ...
    Py = cat(1,Py,py) ;
end
ny = size(Py,1);



% Joint is 6 x 4 ,
% YY YN NY NN for cols
% ai as ah is ih sh rows

% Tbls
% ai as ah is ih sh
% YY YN
% NY NN
%



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
    
    
    
elseif model == 2    % quantum model with 12 parameters used in article
    
    Comp = [1 3];   % arbitrary here, not used
    Inc = [ 2 4];   % arbitrary here, not used
    M2 = ProjQP(2,Comp,Inc,[]);
    
    
    np = 10;    % 10 single 14 both common  16 constrain
    lb = -2*ones(np,1) ;
    ub =  2*ones(np,1) ;
    
    %  use for particle swarm
    fun = @(parm) quant2(parm,M2,Vars,nn,Inc,Py,Fy);
    options1 = optimoptions('particleswarm','SwarmSize',100,'UseParallel',true,'Display','iter','MaxIter',200);
    
    % use for fminunc
    options2 = optimoptions('fminunc','MaxFunEvals',10^10,'MaxIter',1000,'Display','off','Algorithm','quasi-newton');
    
    
    
    PR = zeros(nr,np); ChiR = zeros(nr,1);
    
    for rep = 1:nr
        
        if fit ==1
            [PR(rep,:), ChiR(rep), Px] = particleswarm(fun,np,lb,ub,options1);
            
        elseif fit == 2
            
            % (ai as ah), (is ih), sh
            parm0 = [0.5 -2 0.5 -2 -2 0.35 0.15 0 0 -0.2];
            parm0 = 1*parm0 + .25*randn(1,np);
            % single data set
            [PR(rep,:), ChiR(rep), Px5] = fminunc(@(parm) quant2(parm,M2,[],[],[],Py,Fy),parm0,options2);
            
            [rep ChiR(rep)]
            
        else
            
            if data == 1
                parm =  [0.1924 -1.5400 0.2267 0.2798 -2.0993 0.2769 0.3482 -0.4430 0.3339 -0.5470]; % chi = 42.35
            elseif data == 2
                parm = [0.3163 0.4465 0.2110 -1.6456 -1.9623 0.3717 0.4065 1.1792 0.7627 1.0852 ] ; % chi = 40.45
            elseif data == 3
                parm = [-1.5917 -1.7982 0.5499 -1.8737 -1.5056 -1.7115 0.4467 0.2571 0.4000 -0.0967]; % chi = 37.79
            else 
                % chi = 164.7047
              parm = [-1.7986 -1.7804 0.6462 -1.7966 0.7541 0.7823 0.1393 -0.0200 0.0598 -0.3257 ] ; 
              
            end
            
            
            PR = parm;
            
        end  % pick fit method
        
    end % rep
    
    [Chi, Ind] = min(ChiR);    % pick best fit
    parm = PR(Ind,:);
    
    [Chi,PR, Px] = quant2(parm,M2,Vars,nn,Inc,Py,Fy);

    
    disp('Quantum model')
    disp('Chisquare difference')
    disp(Chi)
    disp('parameters')
    disp(parm)
    
    'predicted top'
    reshape(Py,4,6)' 
                      
    'observed bottom'
    reshape(Px,4,6)'
    % now order is YY NY YN NN
    
end % model

