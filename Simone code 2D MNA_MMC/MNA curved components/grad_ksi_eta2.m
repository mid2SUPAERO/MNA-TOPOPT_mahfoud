%% ***** ksi and eta derivatives with respect to design variables ****** %%
function [dksidx,detadx] = grad_ksi_eta2(x,y,X,n_var)
% ksi & eta are functions of design variables [x0,y0,theta,lx,ly,cr] & (x,y)
% the current point coordinates in global basis
dksidx = zeros(n_var,length(x)) ;
detadx = zeros(n_var,length(x)) ;
x0 = X(1) ;
y0 = X(2) ;
t = X(3) ;
lx = X(4) ;
ly = X(5) ;
cr = X(6) ;
for k = 1:length(x)
    dksi_dx_k = zeros(n_var,1) ;
    deta_dx_k = zeros(n_var,1) ;
    if cr==0
        tau = 10^-6 ;
        s1 = sign(((y(k) - y0 + cos(t)/tau)^2 + (x0 - x(k) + sin(t)/tau)^2)^(1/2) - 1/abs(tau)) ;
        dksi_dx_k(1) = -(2*sign(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0))*cos(t))/lx ;
        dksi_dx_k(2) = -(2*sign(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0))*sin(t))/lx ;
        dksi_dx_k(3) = (2*sign(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0))*(cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0)))/lx ;
        dksi_dx_k(4) = -(2*abs(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0)))/lx^2 ;
        dksi_dx_k(5) = 0 ;
        dksi_dx_k(6) = -(2*s1*(((2*cos(t)*(y(k) - y0 + cos(t)/tau))/tau^2 + (2*sin(t)*(x0 - x(k) + sin(t)/tau))/tau^2)/(2*((y(k) - y0 + cos(t)/tau)^2 + (x0 - x(k) + sin(t)/tau)^2)^(1/2)) - sign(tau)/abs(tau)^2))/ly ;
        angle = pi - 2*atan2(sign(tau)*(1/tau + cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0)), sign(tau)*(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0))) ;
        s2 = sign(angle) ;
        deta_dx_k(1) = (2*sign(cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0))*sin(t))/ly ;
        deta_dx_k(2) = -(2*sign(cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0))*cos(t))/ly ;
        deta_dx_k(3) = -(2*sign(cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0))*(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0)))/ly ;
        deta_dx_k(4) = 0 ;
        deta_dx_k(5) = -(2*abs(cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0)))/ly^2 ;
        deta_dx_k(6) = -(sign(tau)*abs(angle))/(lx*abs(tau)^2) + (2*sign(tau)^2*s2*(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0)))/(tau^2*lx*abs(tau)*(sign(tau)^2*(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0))^2 + sign(tau)^2*(1/tau + cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0))^2)) ;
    else
        s1 = sign(((y(k) - y0 + cos(t)/cr)^2 + (x0 - x(k) + sin(t)/cr)^2)^(1/2) - 1/abs(cr)) ;
        dksi_dx_k(1) = (s1*(2*x0 - 2*x(k) + (2*sin(t))/cr))/(ly*((y(k) - y0 + cos(t)/cr)^2 + (x0 - x(k) + sin(t)/cr)^2)^(1/2)) ;
        dksi_dx_k(2) = -(s1*(2*y(k) - 2*y0 + (2*cos(t))/cr))/(ly*((y(k) - y0 + cos(t)/cr)^2 + (x0 - x(k) + sin(t)/cr)^2)^(1/2)) ;
        dksi_dx_k(3) = (s1*((2*cos(t)*(x0 - x(k) + sin(t)/cr))/cr - (2*sin(t)*(y(k) - y0 + cos(t)/cr))/cr))/(ly*((y(k) - y0 + cos(t)/cr)^2 + (x0 - x(k) + sin(t)/cr)^2)^(1/2)) ;
        dksi_dx_k(4) = 0 ;
        dksi_dx_k(5) = -(2*abs(((y(k) - y0 + cos(t)/cr)^2 + (x0 - x(k) + sin(t)/cr)^2)^(1/2) - 1/abs(cr)))/ly^2 ;
        dksi_dx_k(6) = -(2*s1*(((2*cos(t)*(y(k) - y0 + cos(t)/cr))/cr^2 + (2*sin(t)*(x0 - x(k) + sin(t)/cr))/cr^2)/(2*((y(k) - y0 + cos(t)/cr)^2 + (x0 - x(k) + sin(t)/cr)^2)^(1/2)) - sign(cr)/abs(cr)^2))/ly ;
        angle = pi - 2*atan2(sign(cr)*(1/cr + cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0)), sign(cr)*(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0))) ;
        s2 = sign(angle) ;
        deta_dx_k(1) = -(2*s2*(sin(t)*(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0)) + cos(t)*(1/cr + cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0))))/(lx*abs(cr)*(sign(cr)^2*(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0))^2 + sign(cr)^2*(1/cr + cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0))^2)) ;
        deta_dx_k(2) = (2*s2*(cos(t)*(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0)) - sin(t)*(1/cr + cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0))))/(lx*abs(cr)*(sign(cr)^2*(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0))^2 + sign(cr)^2*(1/cr + cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0))^2)) ;
        deta_dx_k(3) = (2*s2*((cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0))*(1/cr + cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0)) + (cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0))^2))/(lx*abs(cr)*(sign(cr)^2*(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0))^2 + sign(cr)^2*(1/cr + cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0))^2)) ;
        deta_dx_k(4) = -abs(angle)/(lx^2*abs(cr)) ;
        deta_dx_k(5) = 0 ;
        deta_dx_k(6) = -(sign(cr)*abs(angle))/(lx*abs(cr)^2) + (2*sign(cr)^2*s2*(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0)))/(cr^2*lx*abs(cr)*(sign(cr)^2*(cos(t)*(x(k) - x0) + sin(t)*(y(k) - y0))^2 + sign(cr)^2*(1/cr + cos(t)*(y(k) - y0) - sin(t)*(x(k) - x0))^2)) ;
    end
    dksidx(:,k) = dksi_dx_k ;
    detadx(:,k) = deta_dx_k ;
end