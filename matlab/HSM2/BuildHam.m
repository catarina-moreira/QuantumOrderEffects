function [ Ha,Hb ] = BuildHam( P )
% function [ Ha, Hb ] = BuildHam( P )
% Ham transform from |A>*|B> -> |C>

%%%%%%%%%%%%%%%%%%%
% Hamiltonian

ma = P(1);
sa = P(2);
mb = P(3);
sb = P(4);

Ha = [ma  sa ; 
      sa -ma]; 

Hb = [mb  sb   0 ;
      sb  0  sb;
      0   sb  -mb];
  
 
  
end

