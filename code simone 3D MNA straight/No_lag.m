function [FEM_structure] = No_lag(FEM_structure)
DOFs_MAP=FEM_structure.DOFs_MAP;
K_interface=FEM_structure.K_interface;
LOAD_MATRIX=FEM_structure.LOAD_MATRIX;
Recovery_Matrix=FEM_structure.Recovery_Matrix;
U_static=FEM_structure.U_static;
FAN_Center_Recovery_vector=FEM_structure.FAN_Center_Recovery_vector;
U_FAN_static=FEM_structure.U_FAN_static;
ELEMENT=FEM_structure.ELEMENT;
COORD=FEM_structure.COORD;
%% Eliminate Lagrange Multiplier Variables
if sum(DOFs_MAP(:,1)==0)~=0
    nl=sum(DOFs_MAP(:,1)==0);
    Kuu=K_interface(1:end-nl,1:end-nl);
    Kul=K_interface(1:end-nl,end-nl+1:end);
    Kll=K_interface(end-nl+1:end,end-nl+1:end);
    Klu=Kul';
    Fu=LOAD_MATRIX(1:end-nl,:);
    Fl=LOAD_MATRIX(end-nl+1:end,:);
    L0=Kll\Fl;
    Rlu=Kll\Klu;
    K_interface=Kuu-Kul*Rlu;
    LOAD_MATRIX=Fu-Rlu'*Fl;
    Ru=Recovery_Matrix(:,1:end-nl);
    Rl=Recovery_Matrix(:,end-nl+1:end);
    Recovery_Matrix=Ru-Rl*Rlu;
    U_static=U_static+Rl*L0;
    Rfu=FAN_Center_Recovery_vector(:,1:end-nl);
    Rfl=FAN_Center_Recovery_vector(:,end-nl+1:end);
    FAN_Center_Recovery_vector=Rfu-Rfl*Rlu;
    U_FAN_static=U_FAN_static+Rfl*L0;
    DOFs_MAP=DOFs_MAP(DOFs_MAP(:,1)~=0,:);
    %K_interface
    %LOAD_MATRIX
    %Recovery_Matrix
    %FAN_Center_Recovery_vector
    %U_static
    %U_FAN_static
end
%% Elimination of orphan nodes
%analysis
Active_node=unique(ELEMENT(:));
Ic=zeros(size(COORD,1),1);
for k=1:size(COORD,1)
    if any(k==Active_node)
        Ic(k)=find(k==Active_node);
    else
        Ic(k)=0;
    end
end
COORD=COORD(Active_node,:);
ELEMENT=Ic(ELEMENT);
FEM_structure.DOFs_MAP=DOFs_MAP;
FEM_structure.K_interface=K_interface;
FEM_structure.LOAD_MATRIX=LOAD_MATRIX;
FEM_structure.Recovery_Matrix=Recovery_Matrix;
FEM_structure.U_static=U_static;
FEM_structure.FAN_Center_Recovery_vector=FAN_Center_Recovery_vector;
FEM_structure.U_FAN_static=U_FAN_static;
FEM_structure.COORD=COORD;
FEM_structure.ELEMENT=ELEMENT;
