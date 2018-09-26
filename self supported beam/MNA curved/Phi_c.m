%% Forming Phi_i for each component (curved conponents)
function [tmpPhi]=Phi_c(xy,LSgridx,LSgridy,p,fact)
cr=xy(6); % rayon de courbure (inverse de la courbature) 
% changement du rep�re global vers local
theta = xy(3) ;
st = sin(theta) ;
ct = cos(theta) ;
x1 = ct*(LSgridx - xy(1))+st*(LSgridy - xy(2)) ;
y1 = -st*(LSgridx - xy(1))+ct*(LSgridy - xy(2)) ;
% si composant courb�
if ~(cr==0)
    % changement du rep�re local vers centre du cercle
    x2 = x1;
    y2 = y1+1/cr;
    % coordonn�s polaires
    r = sqrt(x2.^2+y2.^2);
    angle = atan2(-sign(cr)*x2,sign(cr)*y2);
    % fix thickness
    tmpPhi = 1-(1/cr*angle).^p/(fact*xy(4))^p-(r-abs(1/cr)).^p/(fact*xy(5))^p ;
% si composant droit
else
    % % variable thickness
    % bb=(xy(5)+xy(4)-2*xy(6))/2/xy(3)^2*x1.^2+(xy(5)-xy(4))/2*x1/xy(3)+xy(6);
    % tmpPhi= -((x1).^p/xy(3)^p+(y1).^p./bb.^p-1);
    % fix thickness
    tmpPhi = 1-(x1).^p/(fact*xy(4))^p-(y1).^p/(fact*xy(5))^p ;
end