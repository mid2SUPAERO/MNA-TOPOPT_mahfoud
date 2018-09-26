%% ********* find the inscribed rectangle inside design domaine ******** %%
% this function finds the inscribed rectangle included in the design domain
% in order to give the initial configuration inside this rectangle, to make
% sure that all components are initially inside the design domain inputs
% are coordinates of nodes and limits of design domaine, this codes take
% account of the shape of the design domain considered here, and it's not
% general at all. outputs are limits of inscribed rectangle edges
function [rect_lim] = inscribed_rectangle(FEM_structure)
coord = FEM_structure.COORD ; % nodes coordinates
n_nod = size(coord,1) ;
limits = FEM_structure.limits ; % limits of design domain
rect_lim(1:2) = limits(1:2) ; 
zmin = limits(5) ; % lower surface z-coordinate
ylist_id = [] ; % list of nodes which belongs to the lower surface 
for i = 1:n_nod
    if coord(i,3)==zmin
        ylist_id = [ylist_id;i] ; %#ok<AGROW>
    end
end
ylist = coord(ylist_id,2) ;
rect_lim(3) = min(ylist) ;
rect_lim(4) = max(ylist) ;
rect_lim(5) = zmin ;
rect_lim(6) = limits(6) ;
end