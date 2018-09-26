%% *************** Forming Phi_i for each component ******************** %%
function [tmpPhi]=TDF_cc(x,meshx,meshy,p)
cr = x(6) ; % courbure
% global to local basis change
st = sin(x(3)) ;
ct = cos(x(3)) ;
x1=ct*(meshx - x(1))+st*(meshy - x(2));
y1=-st*(meshx - x(1))+ct*(meshy - x(2));
if cr == 0 % straight component
    tmpPhi = 1 - (x1).^p/x(4)^p - (y1).^p/x(5)^p ;
else % curved component
    % changement du repère local vers centre du cercle
    x2 = x1;
    y2 = y1+1/cr;
    % coordonnés polaires
    r = sqrt(x2.^2+y2.^2);
    angle = atan2(-sign(cr)*x2,sign(cr)*y2);
    % fix thickness
    tmpPhi = 1-(1/cr*angle).^p/x(4)^p-(r-abs(1/cr)).^p/x(5)^p;
end