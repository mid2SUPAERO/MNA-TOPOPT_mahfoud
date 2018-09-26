%% ********* short script to verify nodes order in an element ********** %%
% in this script we verify how nodes of an element are ordered : if nodes 1
% 2 3 and 4 are coplanar or not
%% nodes of the ith element
i = 20 ; % number of the considered element
n1 = FEM_structure.ELEMENT(i,1) ;
n2 = FEM_structure.ELEMENT(i,2) ;
n3 = FEM_structure.ELEMENT(i,3) ;
n4 = FEM_structure.ELEMENT(i,4) ;
n5 = FEM_structure.ELEMENT(i,5) ;
n6 = FEM_structure.ELEMENT(i,6) ;
n7 = FEM_structure.ELEMENT(i,7) ;
n8 = FEM_structure.ELEMENT(i,8) ;
%% coordinates of nodes
c1 = [FEM_structure.x1(n1);FEM_structure.y1(n1);FEM_structure.z1(n1)] ;
c2 = [FEM_structure.x2(n2);FEM_structure.y2(n2);FEM_structure.z2(n2)] ;
c3 = [FEM_structure.x3(n3);FEM_structure.y3(n3);FEM_structure.z3(n3)] ;
c4 = [FEM_structure.x4(n4);FEM_structure.y4(n4);FEM_structure.z4(n4)] ;
c5 = [FEM_structure.x5(n5);FEM_structure.y5(n5);FEM_structure.z5(n5)] ;
c6 = [FEM_structure.x6(n6);FEM_structure.y6(n6);FEM_structure.z6(n6)] ;
c7 = [FEM_structure.x7(n7);FEM_structure.y7(n7);FEM_structure.z7(n7)] ;
c8 = [FEM_structure.x8(n8);FEM_structure.y8(n8);FEM_structure.z8(n8)] ;
%% vectors definition
v1 = c6-c5 ;
v2 = c7-c5 ;
v3 = c8-c5 ;
%% test of coplanarity
M = [v1 v2 v3] ;
delta = det(M) ;