% Open System Quantum Walk
% Distribution and means at several time points
% compare models PQ, PM , QMc, QMd
% uses continuous time Markov process and generator K



clear
clc
close('all')

Mmethod = 3;  % 1-> T/h t ;  2-> T, t/h ; 3-> T, t
h =   10^-0;   % -2

% gg =   .5;   % (h = 10^-5,  .50  =  5*(10^4)*h  )
gamc = .5;
gamd = gamc;   %  =.5  % wgt on Markov


dim = 21;  % dim = odd no. 21, nt = 100  takes 26 sec
x = (1:dim)';

Max = 5;   % 5
nt = 100; % don't touch
inc = Max/nt;
t = (0:inc:Max);
nt = size(t,2);  % one more including zero


init = (dim+1)./2;
std = 5;   % 5
Psy0 = ( (x-init)./std ).^2;
Psy0 = exp(-Psy0);
Psy0 = Psy0./sqrt(Psy0'*Psy0);


figure
plot(Psy0.^2)
title('Psy0.^2')

% Build Quantum Hamiltonian

a = 0; b = 500/dim; c =0;
pot = a.*x.^2 + b.*x + c;
% a = 1; b = 5;
% pot = a*(x-b).^2;

%  figure
%  plot(x,pot)
% title("potential")

H = zeros(dim,dim);
sig = 35;   %35
sigV = 2*sig*ones(dim-1,1);
H1 = diag(sigV);
H(2:dim,1:(dim-1)) = H1;
H = H + diag(pot);
H = (pi/2)*(H + H')./2;
% H = sparse(H);

% Build Markov Intensity

a = 150; %  100; % off diag top     dec mean pref
b = 200; %  200; % off diag bottom  inc mean pref

gam = -(a+b)*ones(dim,1);
alfa = a*ones(dim,1);    % rate up = (prob up)/h
beta = b*ones(dim,1);    % rate dn = (prob dn)/h

up = alfa(2:dim,1)  ;
dn = beta(1:(dim-1),1);
cn = gam.*ones(dim,1);

A = diag(up);
B = diag(dn);
K = 1*diag(cn);  % 1*
K(1:(dim-1),2:dim) = K(1:(dim-1),2:dim) + A;
K(2:dim,1:(dim-1)) = K(2:dim,1:(dim-1)) + B;
K(1,1) = -K(2,1);
K(dim,dim) = -K(dim-1,dim);
K = sparse(K);

Tdt = expm(K*h);           % transition matrix for small time step
if Mmethod == 1
    G = Tdt/h;
else
    G = Tdt;
end
% [eig(K) eig(G)]

%%%%%%%%%%%%%%%%%%%%  Quantum-Markov walk

rho = Psy0*Psy0';

rhoV = reshape(rho,dim^2,1);

LB = cell(dim,dim);

for j = 1:dim
    for k = 1:dim
        z = 1.*(x == k)*(x==j)' ;
        LB{k,j} = sparse( z );
        % LB{k,j} = sparse( (x == k)*(x==j)' );
    end % k
end  % j

ID = speye(dim);
H = sparse(H);
L1 = -1i*( kron(ID,H) - kron(transpose(H),ID) ) ; % opposite sign of Martinez   

Mc = zeros(dim^2,dim^2); Md = Mc;

for j = 1:dim
    for k = 1:dim
        L = LB{k,j};
        dMc = K(k,j)*( kron(L,L) - .5* (kron(L'*L,ID) + kron(ID,L'*L) )   );
        Mc = dMc + Mc;
        dMd = G(k,j)*( kron(L,L) - .5* (kron(L'*L,ID) + kron(ID,L'*L) )   );
        Md = dMd + Md;
    end % k
end  % j


LMc = (1-gamc)*L1+gamc*Mc;
LMc = sparse(LMc);
LMd = (1-gamd)*L1+gamd*Md;
LMd = sparse(LMd);


%[ Svec, Seig] = eigs(LM,25);


rhoMc = zeros(nt,dim);
rhoMd = rhoMc;
pMark = rhoMc;
pQuan = pMark;

for k = 1:nt
    rhotVc = expm(LMc*t(k))*rhoV;
    rhotc = reshape(rhotVc,dim,dim);
    rhotc = diag(rhotc);
    rhoMc(k,:) = rhotc';
    
    if Mmethod == 2
        rhotVd = expm(LMd*(t(k)/h))*rhoV;
    else
        rhotVd = expm(LMd*t(k))*rhoV;
    end
    % rhotVd = expm(LMd*t(k))*rhoV;
    
    rhotd = reshape(rhotVd,dim,dim);
    rhotd = diag(rhotd);
    rhoMd(k,:) = rhotd';
    pMark(k,:) = (expm(K*t(k))*(Psy0.^2))';
    pQuan(k,:) = (abs(expm(-1i*t(k)*H)*Psy0).^2)';
end



x = x./dim;
meanPc = real(rhoMc*x);
meanPd = real(rhoMd*x);
meanMark = pMark*x;
meanQuan = real(pQuan*x);

figure
plot(t,meanPc,'-r',t,meanPd,'-b',t,meanMark,'-k',t,meanQuan,'-g','linewidth',2);
xlabel('Time')
ylabel('Mean Evidence')
axis( [ 0 Max .25 1] )
legc = strcat('QMc ','- ',num2str(gamc));
legd = strcat('QMd ','- ',num2str(gamd));
legend(legc,legd,'Mark','Quant')
title( " Mean evidence plotted against time")


[r,v] = min(real(rhoMc));
disp('r is min rhoMc')
disp(r)

disp('      Pc      Pd       Markov     Quan')
disp(real([meanPc meanPd   meanMark  meanQuan]))

figure
tt = nt;
sz = size(tt,2);
plot(x*ones(1,sz),real(rhoMc(tt,:))','-r', x*ones(1,sz),real(rhoMd(tt,:))','-b', x,pMark(tt,:),'-k','linewidth',2)
xlabel('Evidence value')
ylabel('Probability')
legend('Quantum C','Quantum D','Markov')
title('Probability distribution across states at last time point')

