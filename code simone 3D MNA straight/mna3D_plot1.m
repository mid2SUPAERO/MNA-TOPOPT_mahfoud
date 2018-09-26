%% ********************* plot MNA components in 3D ********************* %%
% this function draws component plots under MNA framework in 3D, inputs are
% design variables, the Topology Description function is used to plot the
% the contour of the structure in 3D
function mna3D_plot1(formulation,FEM_structure)
x = formulation.x ; % design vector
n_var = FEM_structure.n_var ; % total number of variables
n_var_c = FEM_structure.n_var_c ; % number of variables per components
n_c = n_var/n_var_c ; % number of components
