%% *************** Forming Phi_i for each component ******************** %%
function [tmpPhi]=TDF(x,meshx,meshy,p)
st = sin(x(3)) ;
ct = cos(x(3)) ;
x1=ct*(meshx - x(1))+st*(meshy - x(2));
y1=-st*(meshx - x(1))+ct*(meshy - x(2));
tmpPhi= 1 - (x1).^p/x(4)^p - (y1).^p/x(5)^p ;
end