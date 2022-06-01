function M = ProjQP(nv,Comp,Inc,Hdim)
% function M = ProjQP(nv,Comp,Inc)

% create projectors for quantum model


na = nv(1);  % no. of values for self rating
nb = nv(2);  % no. of values for other rating
nc = nv(3);  % 1 
nd = nv(4);  % 1

Sdima = Hdim/na;   % 11
Sdimb = Hdim/nb;   % 11

cq = cell(1,2);


ca = cell(1,na);  % 9 projectors
for j=1:na
    ca{j} = (1:Sdima) + ((j-1)*Sdima);
end

cb = cell(1,nb);
for j=1:nb
    cb{j} = (1:Sdimb) + ((j-1)*Sdimb);
end

cq{1} = ca;
cq{2} = cb;

MA = BuildProj(cq{Inc(1)});
MB = BuildProj(cq{Inc(2)});
M = {MA,MB,1,1};

