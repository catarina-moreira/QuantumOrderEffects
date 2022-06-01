function M = ProjQP(nv,~,~,~)
%function M = ProjQP(nv,Comp,Inc,Hdim)
% create projectors for quantum model


ca = cell(1,nv);
for j=1:nv
    ca{j} = j;
end

M = cell(1,2);
[Mx,My] = BuildComp2(ca,ca);
M([1 2]) = {Mx, My};  


% YY
% YN
% NY
% NN


% |x*y>
% Mx(1)  picks yes to first cubit attibute  YY+YN
% My(1)  picks yes to second  YY + NY
