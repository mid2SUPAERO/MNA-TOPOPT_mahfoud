%% *********** global cartesean to local polar basis change ************ %%
function [ksi,eta] = changebasis3(X,Y,dsgn)
x0 = dsgn(1) ;
y0 = dsgn(2) ;
theta = dsgn(3) ;
lx = dsgn(4) ;
ly = dsgn(5) ;
cr = dsgn(6) ;
x = X(:) ;
y = Y(:) ;
% global to local cartesean basis change
st=sin(theta);
ct=cos(theta);
x1=ct*(x - x0)+st*(y - y0);
y1=-st*(x - x0)+ct*(y - y0);
% si composant courbé
if ~(cr==0)
    % changement du repère local vers centre du cercle
    x2=sign(cr)*x1;
    y2=sign(cr)*(y1+1/cr);
    % coordonnés polaires
    r=sqrt(x2.^2+y2.^2);
    phi=atan2(y2,x2);
    % normalization
    ksi = 2/ly*abs(r-1/abs(cr)) ;
    eta = 1/lx*abs((pi-2*phi)/cr) ;
% si composant droit
else
    ksi = 2/lx*abs(x1) ;
    eta =2/ly*abs(y1) ;
end