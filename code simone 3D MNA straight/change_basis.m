%% ******************* local to global basis change ******************** %%
%this function transforms coordinates of a point from local to global basis
function [x,y,z] = change_basis(ksi,eta,mu,x,y,z,X)
% angles and lengths short-hand
c1 = cos(X(7)) ; c2 = cos(X(8)) ; c3 = cos(X(9)) ;
s1 = sin(X(7)) ; s2 = sin(X(8)) ; s3 = sin(X(9)) ;
l1 = X(4) ; l2 = X(5) ; l3 = X(6) ;
% rotation matrix 
R = [     c2*c3          -c2*s3          s2   ;...
     s1*s2*c3+c1*s3  -s1*s2*s3+c1*c3   -s1*c2 ;...
    -c1*s2*c3+s1*s3   c1*s2*s3+s1*c3    c1*c2 ] ;
% global coordinates
V_global = R\[ksi*l1;eta*l2;mu*l3]+[x;y;z] ;
x = V_global(1,:) ; y = V_global(2,:) ; z = V_global(3,:) ;
