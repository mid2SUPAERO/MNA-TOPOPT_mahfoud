function [FEM_structure]=mesh_refinement_3D(FEM_structure,nfin)
% make nfin refinement on mesh breaking each element in 4
ELEMENT=FEM_structure.ELEMENT;
COORD=FEM_structure.COORD;
FACES=FEM_structure.FACES;
NODE_SET=FEM_structure.NODE_SET;
FEM_structure.Old_coord=COORD;
FEM_structure.Old_element=ELEMENT;
%% refinement
ELEMENTc=ELEMENT;
Element_id=1:size(ELEMENTc,1);
Pu=cell(nfin-1,1);
for l = 1:nfin
    %     [Pu{l,1}] = Prepcoarse(nely/2^(l-1),nelx/2^(l-1));
    [COORD,ELEMENT,NODE_SET,FACES,Pu{nfin-l+1,1}]=mesh_refinement3Dpar(COORD,ELEMENT,NODE_SET,FACES);
    Element_id=reshape(repmat(Element_id,8,1),[],1).';
end
FEM_structure.COORD=COORD;
FEM_structure.ELEMENT=ELEMENT;
FEM_structure.NODE_SET=NODE_SET;
FEM_structure.Element_id=Element_id;
FEM_structure.Pu=Pu;
FEM_structure.FACES=FACES;