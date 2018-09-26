%% ******************* create initial configuration ******************** %%
% this function initialize the design configuration, from number of
% components desired in each direction and dimensions of design domaine,
% and gives a regular distribution of components in the design domaine
function [design_var] = initial_design(FEM_structure)
n_var_c = FEM_structure.n_var_c ;
nX = FEM_structure.nX ;
nY = FEM_structure.nY ;
nZ = FEM_structure.nZ ;
design_var = zeros(n_var_c,nX*nY*nZ) ;
l1 = FEM_structure.limits(2)-FEM_structure.limits(1) ;
l2 = FEM_structure.limits(4)-FEM_structure.limits(3) ;
l3 = FEM_structure.limits(6)-FEM_structure.limits(5) ;
d = ((0.9*l1*l2*l3)/(nX*nY*nZ))^(1/3) ;
design_var(4,:) = d ;
design_var(5,:) = d ;
design_var(6,:) = d ;
[rect_lim] = inscribed_rectangle(FEM_structure) ;
% x = linspace(FEM_structure.limits(1),FEM_structure.limits(2),nX+2) ;
% x = setdiff(x,[x(1) x(nX+2)]) ;
% y = linspace(FEM_structure.limits(3),FEM_structure.limits(4),nY+2) ;
% y = setdiff(y,[y(1) y(nY+2)]) ;
% z = linspace(FEM_structure.limits(5),FEM_structure.limits(6),nZ+2) ;
% z = setdiff(z,[z(1),z(nZ+2)]) ;
x = linspace(rect_lim(1),rect_lim(2),nX+2) ;
x = setdiff(x,[x(1) x(nX+2)]) ;
y = linspace(rect_lim(3),rect_lim(4),nY+2) ;
y = setdiff(y,[y(1) y(nY+2)]) ;
z = linspace(rect_lim(5),rect_lim(6),nZ+2) ;
z = setdiff(z,[z(1),z(nZ+2)]) ;
[X,Y,Z] = meshgrid(x,y,z) ;
design_var(1,:) = X(:) ;
design_var(2,:) = Y(:) ;
design_var(3,:) = Z(:) ;
design_var = reshape(design_var,FEM_structure.n_var,1) ;
end