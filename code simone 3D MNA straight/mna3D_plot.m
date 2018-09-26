%% ********************* plot MNA components in 3D ********************* %%
% this function draws component plots under MNA framework in 3D, inputs are
% design variables, the Matlab function "patch" is used to plot the
% polygons in 3D
function mna3D_plot(formulation,FEM_structure)
x = formulation.x ; % design vector
n_var = FEM_structure.n_var ; % total number of variables
n_var_c = FEM_structure.n_var_c ; % number of variables per components
n_c = n_var/n_var_c ; % number of components
% local coordinates of components vertexes
scale = 1 ; % scaling factor to have good plots
ksi = scale*[-1 1 1 -1 -1 1 1 -1] ; 
eta = scale*[-1 -1 1 1 -1 -1 1 1] ;
mu = scale*[-1 -1 -1 -1 1 1 1 1] ;
figure(2) ;
clf
for i = 1:n_c
    j = n_var_c*(i-1)+1 ;
    X = repmat(x(j),1,8) ; Y = repmat(x(j+1),1,8) ; Z = repmat(x(j+2),1,8);
    [X,Y,Z] = change_basis(ksi,eta,mu,X,Y,Z,x(j:j+8)) ;
    % faces of component
    face1 = 1:4 ; face2 = 5:8 ; face3 = [1 2 6 5] ;
    face4 = [2 3 7 6] ; face5 = [3 4 8 7] ; face6 = [4 1 5 8] ;
    patch(X(face1),Y(face1),Z(face1),[0 0 1]) ;
    patch(X(face2),Y(face2),Z(face2),[0 0 1]) ;
    patch(X(face3),Y(face3),Z(face3),[0 0 1]) ;
    patch(X(face4),Y(face4),Z(face4),[0 0 1]) ;
    patch(X(face5),Y(face5),Z(face5),[0 0 1]) ;
    patch(X(face6),Y(face6),Z(face6),[0 0 1]) ;
end
