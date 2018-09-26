%% ********* gradient of ksi, eta and mu wrt design variables ********** %%
% the purpose of this script is to get analytical expressions of gradients
% of ksi, eta and mu with respect to design variables. Maltlab Symbolic
% math toolbox is used. the expressions calculated are needed for density
% gradient calculation.
clear ; close all ;
%% variables definition
syms x y z ; % coordinates of current point in global basis
syms x0 y0 z0 ; % coordinates of mass node (design variables)
syms l1 l2 l3 ; % half lenghts of component (design variables)
syms t1 t2 t3 ; % orientation angles (design variables)
syms ksi eta mu ; % coordiantes of current point in local basis
%% global to local basis change
% sines and cosines short-hand
c1 = cos(t1) ; c2 = cos(t2) ; c3 = cos(t3) ;
s1 = sin(t1) ; s2 = sin(t2) ; s3 = sin(t3) ;
% rotation matrix 
R = [     c2*c3          -c2*s3          s2   ;...
     s1*s2*c3+c1*s3  -s1*s2*s3+c1*c3   -s1*c2 ;...
    -c1*s2*c3+s1*s3   c1*s2*s3+s1*c3    c1*c2 ] ;
% local coordiantes
V_local = R*[x-x0;y-y0;z-z0] ;
ksi = abs(V_local(1))/l1 ; eta = abs(V_local(2))/l2 ;mu=abs(V_local(3))/l3;
%% gradient calculation
grad_ksi = gradient(ksi,[x0 y0 z0 l1 l2 l3 t1 t2 t3]) ;
grad_eta = gradient(eta,[x0 y0 z0 l1 l2 l3 t1 t2 t3]) ;
grad_mu = gradient(mu,[x0 y0 z0 l1 l2 l3 t1 t2 t3]) ;