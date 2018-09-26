%% ******************* global to local basis change ******************** %%
% this function changes the coordinates of element centers from global to
% local basis in order to compute density field and its gradient
function [ksi,eta,mu] = basis_change(xc,yc,zc,X)
% sines and cosines short-hand
c1 = cos(X(7)) ; c2 = cos(X(8)) ; c3 = cos(X(9)) ;
s1 = sin(X(7)) ; s2 = sin(X(8)) ; s3 = sin(X(9)) ;
% rotation matrix 
R = [     c2*c3          -c2*s3          s2   ;...
     s1*s2*c3+c1*s3  -s1*s2*s3+c1*c3   -s1*c2 ;...
    -c1*s2*c3+s1*s3   c1*s2*s3+s1*c3    c1*c2 ] ;
% local coordiantes
V_local = R*[xc-X(1);yc-X(2);zc-X(3)] ;
l1 = X(4) ; l2 = X(5) ; l3 = X(6) ;
ksi = abs(V_local(1))/l1 ; eta = abs(V_local(2))/l2 ;mu=abs(V_local(3))/l3;