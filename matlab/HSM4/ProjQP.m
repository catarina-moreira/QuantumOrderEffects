function M = ProjQP(nv,Comp,Inc,~)
% function M = ProjQP(nv,Comp,Inc,Hdim)
% create projectors for quantum model

% Create Projectors
na = nv(1);
nb = nv(2);
nc = nv(3);
nd = nv(4);


cq = cell(1,3);

ca = cell(1,na);
for j=1:na
    ca{j} = j;
end

cb = cell(1,nb);
for j=1:nb
    cb{j} = j;
end


cq{1} = ca;
cq{2} = cb;

[MA, MB] = BuildComp2(cq{Comp(1)},cq{Comp(2)});
MC = MB;
M = {MA,MB,MC,1};

