function [ H ] = BuildHam( P,pos )
% function [ H ] = BuildHam( P,pos )
% Ham transform  |A>*|B> 
% pos refers first or sec cubit in 4-D space
% Ham transform  |A>*|B>   A is first pos, B is sec
%%%%%%%%%%%%%%%%%%%
% Ham


%  pos     first att           sec att
mu = (pos==1)*HTan(P(2)) + (pos==2)*HTan(P(3));      
gam = P(1); 

Ha =(1/sqrt(1+(mu.^2)))*[mu 1; 1 -mu] ;

% rotate for first attribute

H1 = kron(Ha,eye(2)) ;

% rotate for second attribute

H2 = kron(eye(2),Ha) ;


Hc = (1/sqrt(2))*( kron([1 1;1 -1],[1 0; 0 0])...
    + kron([-1 1;1 1],[0 0; 0 1]) );
Hc = -gam*Hc;

% pick j=1 first att, 2 = second


H = (pos==1)*(H1) + (pos==2)*(H2+Hc);

% H = H + Hc;


end

