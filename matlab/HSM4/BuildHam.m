function H = BuildHam( P )
% function H = BuildHam( P )
% Ham transform from |A>*|B> -> |C>

%%%%%%%%%%%%%%%%%%%
% Hamiltonian

ma = P(1);
mb = P(2);
s = P(3);
gam = P(4);

H = zeros(4,4);

Ha = [ma  s ; 
      s -ma]; 

Hb = [mb  s   ;
      s -mb];
  
H(1:2,1:2) = Ha;
H(3:4,3:4) = Hb;


Hc = [ 1  0  1  0 ;
       0 -1  0  1 ; 
       1  0 -1  0 ;
       0  1  0  1];

H = H - (gam/sqrt(2))*Hc;
  
end

