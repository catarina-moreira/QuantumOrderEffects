function [ TAB ] = TwoWayM( A,B,Psy,na,nb)
% function [ Tbl ] = TwoWayM( MA,MB, Psy,na,nb )
% create two way contingency table from Q model

sz = size(Psy,1);
J = ones(sz,1);

TAB = zeros(na,nb);
for j = 1:na
    for k = 1:nb
        M1 = squeeze(A(j,:,:));
        M2 = squeeze(B(k,:,:));
        TAB(j,k) = J'*M2*M1*Psy;
    end
end

end

