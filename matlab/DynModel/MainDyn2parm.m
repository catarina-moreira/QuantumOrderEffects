% Main HSM Dynamic model program
% two parameter models 5 categories

clear
clc

%%

% choose data to fit
% Choose DynConf
data = 1;   % this is the only one for this example

% choose model to fit
model = 1;   %  1= Markov , 2 =quant

% choose to fit data or to check fit
fit = 0 ;  % 0 = fit off compute predictions only, 1 = use fminun
% number of replications of starting points (especially needed for fminunc)
reps = (fit==1).*5 + (fit==0).*1;

% select subj
    subj = 2;      % select subj 

% Quantum dyn model for rating scales
% In this example
% Fit Confidence Ratings at time 1 and then at time 2
% vary the time 1 , time 2, points across 3 conditions

% your input needed
% ******************************
%  enter your # of values for each variable
ns=5;  % input # values of time 1 ratings
no=ns;  % input # values of time 2 ratings 
% (same as time 1 because we are using same measurement twice)
nc=1;  % input #   fix to one
nd=1;  % input #   fix to one
nv = [ns no nc nd];
Jdim = prod(nv);
Hdim = 100; % 100 states used in Hilbert space 
% this can be changed, but Hdim must be > ns

%%

% your input needed
% ******************
% read in your raw frequency tables into a cell called Tbls


if data == 1
    % you need to input your raw freq tables
    
    % this example has a cell with three joint freq tables
    % prob transit row at t=.5 to col at t=1.5
    % prob transit row at t=1.5 to col at t=2.5
    % prob transit row at t=.5 to col at t=2.5
    % load DynConfData
    
    load DblConfDat5
    
    % this file contains a cell(11,3)  11 subjects , 3 tables
    % ( Dbl conf data, summed over coherence )
    % rows are subjects, columns are tables
    % JF1 %  FS{j,1}  row t1 col t2 for subj j
    % JF2 %  FS{j,2}  row t2 col t3 for subj j
    % JF3 %  FS{j,3}  row t1 col t3 for subj j
    
    Tbls = FS(subj,:);
    % 1 x 3 cell for subj j
    % contain 3 tables  {JF1, JF2, JF3}
    
else   % include your data here
    
end   % type of data


% This defines the variables in the tables
% This defines three tables for the JF1, JF2, JF3 
% 1 = first conf measure , 2 = second conf measure
Vars = { [1 2], [1 2], [1 2] };

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

%%

if model == 1

    % Markov
    Comp = 1;
    Inc = 2;
    M2 = ProjQP(nv,Comp,Inc,Hdim);  % creates projectors for rating levels
    
    np = 2;
    % options = optimset('Display','iter','MaxFunEvals',2000,'MaxIter',1000,'TolX',1e-6);
    options = optimset('Display','off','MaxFunEvals',2000,'MaxIter',1000,'TolX',1e-3);
    
    
    ChiV = zeros(reps,1);
    ParmM = zeros(reps,np);
    
    std = 1;    % amount of jitter
    
    for n = 1:reps
        
        disp(n)
        jitter = std*randn(1,np);
        %         pot  diff  start time intervals
        parm0 =  [ 200 200   ];
        parm0 = parm0 + jitter;      
          
        if fit == 1
            
            [parm,Chi] = fminsearch(@(parm) Markov(parm,M2,Vars,nn,Inc,Py,Fy), parm0, options);
            [n Chi/1000]
            ChiV(n) =  Chi;
            ParmM(n,:) =  parm;
        else
            
            
%             % subj 1 , ns = 5 levels  
%             parm = [998.0560 1.0273e+03  ];  % wd=3
%             ChiV =  2.1326e+03;
%             parm = [  987.3    1016.4] ;   % 5 equal space 
%             ChiV = 2.1293e+03;
%             parm = [377.5646  389.0362 ] ;
%             ChiV = 921.9511
              
%             % subj 2 
%             parm = [901.1074  956.7971 ] ;  % 3
%             ChiV =  926.7278; 
%             parm = [892.5386  947.9177 ] ;  % 5
%             ChiV =  925.7225; 
              parm = [239.1234  258.1181 ] ;
%             ChiV = 482.6808;

%             % subj 3
%             parm = [801.2498  827.3024 ] ; % 3
%             ChiV =  1.5442e+03 ;
%             parm = [790.2075  816.0430 ] ; % 5
%             ChiV =  1.5416e+03 ;

%             % subj 4
%             parm = [ 2.6210    2.9448 ] ;  % 3
%             ChiV = 643.5592 ;
%             parm = [ 0.7733    1.0553 ] ;  % 5
%             ChiV = 619.3011 ;

%             % subj 5
%             parm = [1.0961e+03 1.1394e+03] ;
%             ChiV = 1.1590e+03;
%             parm = [1.0862    1.1293] ;
%             ChiV = 1.1577e+03;

%             % subj 6
%             parm = [803.5827  804.9238 ] ;
%             ChiV =  2.1471e+03;
%             parm = [ 793.0758  794.3997 ] ;
%             ChiV =  2.1448e+03;

%             % subj 7
%             parm = [1249.3    1295.7 ] ;
%             ChiV =  3.3272e+03;
%             parm = [1239.7    1285.9 ] ;
%             ChiV =  3.3250e+03;

%             % subj 8 
%             parm = [1071.8    1097.1];
%             ChiV = 2.7183e+03;
%             parm = [1060.5    1085.6];
%             ChiV = 2.7148e+03;

%             % subj 9 
%             parm = [918.8799  951.9521 ] ;
%             ChiV = 3.0533e+03;
%             parm = [908.1595  941.0087 ] ;
%             ChiV = 3.0495e+03;

%             % subj 10
%             parm = [1460.8    1488.9 ] ;
%             ChiV = 2.4614e+03 ;
%             parm = [1451.2    1479.2 ] ;
%             ChiV = 2.4607e+03 ;

            % subj 11
%             parm = [308.0507  314.3267 ] ;
%             ChiV = 2.3577e+03;
%              parm = [300.3970  306.6058 ] ;
%              ChiV = 2.3493e+03;
            
            ParmM = parm;
        end  % fit method
    end  % reps
    
    [Chi, Ind] = min(ChiV);    % pick best fit Index
    parm = ParmM(Ind,:);
    
    [Chi,parm, Px] = Markov(parm,M2,Vars,nn,Inc,Py,Fy);
    
    disp('Markov model')
    disp('Chisquare difference')
    disp(Chi)
    
    % save MFitResults parm PM Fdata
    
    disp('Parms')
    disp(parm)
    
    
    
          'Predicted'
         Px = reshape(Px,ns,ns,3)  % prediction
%          squeeze(sum(Px,1))
          'observed'
          Py = reshape(Py,ns,ns,3)  % observed
%        % squeeze(sum(Py,1))
    
    

%% ******************************** 

else
    
    Comp = 1;
    Inc = 2;
    M2 = ProjQP(nv,Comp,Inc,Hdim);
    
    np = 2;
    
    %  use for particle swarm
    fun = @(parm) quant2(parm,M2,Vars,nn,Inc,Py,Fy);
    options1 = optimoptions('particleswarm','SwarmSize',100,'UseParallel',true,'Display','iter','MaxIter',200);
    lb =  0*ones(np,1) ;  % only for swarm
    ub =  500*ones(np,1) ;
    
    % options = optimset('Display','iter','MaxFunEvals',2000,'MaxIter',1000,'TolX',1e-6);
    options = optimset('Display','off','MaxFunEvals',2000,'MaxIter',1000,'TolX',1e-3);
    
    
    ChiV = zeros(reps,1);
    ParmM = zeros(reps,np);
    
    std = 50;    % amount of jitter
    
    for n = 1:reps
        
        disp(n)
        jitter = std*randn(1,np);
        %         pot  diff  start time intervals
        parm0 =  [    103.8009  227.9163  ];
        
        parm0 = parm0 + jitter;
        
        if fit == 1
            
             [parm,Chi] = particleswarm(fun,np,lb,ub,options1);
           % [parm,Chi] = fminsearch(@(parm) quant2(parm,M2,Vars,nn,Inc,Py,Fy), parm0, options);
            [n Chi/1000]
            ChiV(n) =  Chi;
            ParmM(n,:) =  parm;
            
        else
            
            
            % subj 1 , ns = 5 levels wd  = 3
            
%                parm = [100.0669  208.7825 ] ;    % equal space
%                ChiV = 1.4422e+03;
%                parm = [  103.8009  227.9163 ];  % unequal space
%                ChiV = 396.0744
              
            % subj 2, ns = 5 levels
%             parm = [ 162.8878  134.8060 ];
%             ChiV = 928.0008;
%             parm = [164.1413  140.1246 ] ;
%             ChiV = 363.5141;

%             subj 3 ns = 5
%             parm = [96.9594  210.4032];
%             ChiV = 1.0849e+03;

%             subj 4 
%               parm = [ 10.8725    2.3877] ;
%               ChiV = 669.2089 ;
%               
%         %      subj 5
%               parm = [ 99.2861  209.1001 ] ;
%               ChiV = 981.3916 ;
%               
%               subj 6
%                 parm = [ 99.3060  282.2715] ;
%                 ChiV = 1.2936e+03 ;

%               subj 7
%                 parm = [99.4297  209.0913 ] ;
%                 ChiV = 2.7297e+03 ;

%               subj 8
%                 parm = [99.8652  208.7859 ] ;
%                 ChiV =  1.9008e+03;
                
%               subj 9 
%                 parm = [99.2779  209.1247 ] ;
%                 ChiV = 2.4002e+03;

%               subj 10
%                 parm = [62.5422  288.4035];
%                 ChiV = 1.9773e+03;
                
%               subj 11
%                 parm =   [  99.2792  281.3576];
%                 ChiV = 1.8635e+03;  
%                 
                
            ParmM = parm;
            
        end  % fit method
    end  % reps
    
    [Chi, Ind] = min(ChiV);    % pick best fit Index
    parm = ParmM(Ind,:);
    
    [Chi,parm, Px] = quant2(parm,M2,Vars,nn,Inc,Py,Fy);
    
    disp('Quantum model')
    disp('Chisquare difference')
    disp(Chi)
    
    % save QFitResults parm PM Fdata
    
    disp('Parms')
    disp(parm)
  
    
    'Predicted'
    Px = reshape(Px,ns,ns,3)  % prediction
    %  squeeze(sum(Px,2))
    'observed'
    Py = reshape(Py,ns,ns,3)  % observed
    %  squeeze(sum(Py,1))
    
    
end  % type of model fit
