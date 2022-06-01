function M = ProjQP(nv,~,~,Hdim)
% function M = ProjQP(nv,Comp,Inc,Hdim)
% create projectors for quantum model

scaling = 2;

% if na=5 then build 5 projectors
na = nv(1);  % no. of values for first rating, e.g. 5

if scaling == 1
    % equally spaced states per rating 
    Sdima = Hdim/na;
    Sdima = ones(5,1)*Sdima;
else
    % unequally spaced states for ratings 
    % set spacing here
    Sdima = [40 5 10 5 40]; 
end


ca = cell(1,na);
for j=1:na
    ca{j} = (1:Sdima(j)) + sum(Sdima(1:(j-1))  );
end



MA = BuildProj(ca);
M = {MA,MA};

