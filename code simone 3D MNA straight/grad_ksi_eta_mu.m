%% ********* gradient of ksi, eta and mu wrt design variables ********** %%
% this function computes the analytical gradients of local coordiantes
% ksi eta and mu with respect to design variables. the symbolic expressions
% are token from the code "ksi_eta_mu_grad.m" and copied in this function,
% the gradients are calculated in different points so that grad_ksi,
% respectively grad_eta, grad_mu is a length(ksi) by number of considered
% design variables matrix. input variable neigh_g is a matrix whose rows 
% are xc yc and zc global coordinates of neighborhood element centers. X is
% a vector containing design variables of a component.  
function [grad_ksi,grad_eta,grad_mu] = grad_ksi_eta_mu(neigh_g,X)
ln = length(neigh_g(1,:)) ;
x = neigh_g(1,:) ; y = neigh_g(2,:) ; z = neigh_g(3,:) ;
grad_ksi = zeros(ln,9) ;
grad_eta = zeros(ln,9) ;
grad_mu = zeros(ln,9) ;
x0 = X(1); y0 = X(2); z0 = X(3); l1 = X(4); l2 = X(5); l3 = X(6);
t1 = X(7) ; t2 = X(8) ; t3 = X(9) ;
A = sin(t2)*(z - z0) + cos(t2)*cos(t3)*(x - x0) - cos(t2)*sin(t3)*(y - y0);
grad_ksi(:,1) = -(sign(A)*cos(t2)*cos(t3))/l1 ;
grad_ksi(:,2) = (sign(A)*cos(t2)*sin(t3))/l1 ;
grad_ksi(:,3) = -(sign(A)*sin(t2))/l1 ;
grad_ksi(:,4) = -abs(A)/l1^2 ;
grad_ksi(:,5) = 0 ;
grad_ksi(:,6) = 0 ;
grad_ksi(:,7) = 0 ;
grad_ksi(:,8) = (sign(A).*(cos(t2)*(z - z0) - cos(t3)*sin(t2)*(x - x0) + ...
    sin(t2)*sin(t3)*(y - y0)))/l1 ;
grad_ksi(:,9) = -(sign(A).*(cos(t2)*cos(t3)*(y - y0) + ...
    cos(t2)*sin(t3)*(x - x0)))/l1 ;
B = (cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2))*(x - x0) + ...
    (cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3))*(y - y0)...
    - cos(t2)*sin(t1)*(z - z0) ;
grad_eta(:,1) = -(sign(B)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2)))/l2 ;
grad_eta(:,2) = -(sign(B)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)))/l2 ;
grad_eta(:,3) = (sign(B)*cos(t2)*sin(t1))/l2 ;
grad_eta(:,4) = 0 ;
grad_eta(:,5) = -abs(B)/l2^2 ;
grad_eta(:,6) = 0 ;
grad_eta(:,7) = -(sign(B).*((sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2))*...
    (x - x0) + (cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3))*(y - y0)...
    + cos(t1)*cos(t2)*(z - z0)))/l2 ;
grad_eta(:,8) = (sign(B).*(sin(t1)*sin(t2)*(z - z0) + cos(t2)*cos(t3)...
    *sin(t1)*(x - x0) - cos(t2)*sin(t1)*sin(t3)*(y - y0)))/l2 ;
grad_eta(:,9) = (sign(B).*((cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3))*...
    (x - x0) - (cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2))*(y - y0)))/l2 ;
C = (sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2))*(x - x0) + (cos(t3)*...
    sin(t1) + cos(t1)*sin(t2)*sin(t3))*(y - y0) + cos(t1)*cos(t2)*(z-z0) ;
grad_mu(:,1) = -(sign(C)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2)))/l3 ;
grad_mu(:,2) = -(sign(C)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)))/l3 ;
grad_mu(:,3) = -(sign(C)*cos(t1)*cos(t2))/l3 ;
grad_mu(:,4) = 0 ;
grad_mu(:,5) = 0 ;
grad_mu(:,6) = -abs(C)/l3^2 ;
grad_mu(:,7) = (sign(C).*((cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2))*...
    (x - x0) + (cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3))*(y - y0)...
    - cos(t2)*sin(t1)*(z - z0)))/l3 ;
grad_mu(:,8) = -(sign(C).*(cos(t1)*sin(t2)*(z - z0) + cos(t1)*cos(t2)*...
    cos(t3)*(x - x0) - cos(t1)*cos(t2)*sin(t3)*(y - y0)))/l3 ;
grad_mu(:,9) = (sign(C).*((cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3))*...
    (x - x0) - (sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2))*(y - y0)))/l3 ;
end