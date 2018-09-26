clear
close all
load('substructure_results')
load('lin_pert10.dat.mat')
X_dis=unique(intercoord(:,2));
nx=length(X_dis);
dz=mean(diff(X_dis));
nz=round(nx/2);
Lz=nz*dz;
%the point at the same x have the same 1 and 3 diplacement
%projection_matrix
Pr=zeros(size(K_interface,1),2*nx+sum(DOFs_MAP(:,1)==0));
new_DOFs_MAP=zeros(2*nx+sum(DOFs_MAP(:,1)==0),4);
new_DOFs_MAP(:,4)=(1:2*nx+sum(DOFs_MAP(:,1)==0)).';
for k=1:nx
    x=X_dis(k);
    concerned_nodes=intercoord(intercoord(:,2)==x,1);
    for dof=1:2
        DOF=[1;3];
        [Ind]=ismember(DOFs_MAP(:,1),concerned_nodes);
        concerned_dofs=DOFs_MAP(Ind,4);
        %choose a retained y=0;
        toll=1e-9;
        central_nod=intercoord(intercoord(:,2)==x&abs(intercoord(:,3))<=toll,1);
        new_DOFs_MAP((k-1)*2+dof,[1,3])=DOFs_MAP(DOFs_MAP(:,1)==central_nod&DOFs_MAP(:,3)==DOF(dof),[1,3]);
        new_DOFs_MAP((k-1)*2+dof,2)=k;
        Pr(Ind,(k-1)*2+dof)=1;
    end
end
%add DOF for LAGRANGIAN multipliers
Pr(end-sum(DOFs_MAP(:,1)==0)+1:end,end-sum(DOFs_MAP(:,1)==0)+1:end)=eye(sum(DOFs_MAP(:,1)==0));
new_DOFs_MAP(end-sum(DOFs_MAP(:,1)==0)+1:end,1:3)=DOFs_MAP(DOFs_MAP(:,1)==0,1:3);
K_int_new=Pr.'*K_interface*Pr;
LOAD_MATRIX_new=Pr.'*LOAD_MATRIX;
Recovery_Matrix_new=Recovery_Matrix*Pr;
FAN_Center_Recovery_vector_new=FAN_Center_Recovery_vector*Pr;
%
Z_0=intercoord(intercoord(:,1)==central_nod,4);
Z_dis=Z_0:dz:Z_0+Lz;
nz=length(Z_dis);
new_DOFs_MAP(end-sum(DOFs_MAP(:,1)==0)+1:end,2)=0;
new_DOFs_MAP(end-sum(DOFs_MAP(:,1)==0)+1:end,4)=(1:sum(DOFs_MAP(:,1)==0)).'+2*nz*nx;
Kcc=K_int_new(1:end-sum(DOFs_MAP(:,1)==0),1:end-sum(DOFs_MAP(:,1)==0));
Klc=K_int_new(end-sum(DOFs_MAP(:,1)==0)+1:end,1:end-sum(DOFs_MAP(:,1)==0));
Kll=K_int_new(end-sum(DOFs_MAP(:,1)==0)+1:end,end-sum(DOFs_MAP(:,1)==0)+1:end);
K_def=zeros(2*nx*nz+sum(DOFs_MAP(:,1)==0));
K_def(1:2*nx,1:2*nx)=Kcc;
K_def(1:2*nx,end-sum(DOFs_MAP(:,1)==0)+1:end)=Klc.';
K_def(end-sum(DOFs_MAP(:,1)==0)+1:end,1:2*nx)=Klc;
K_def(end-sum(DOFs_MAP(:,1)==0)+1:end,end-sum(DOFs_MAP(:,1)==0)+1:end)=Kll.';
K_def=sparse(K_def);
LOAD_MATRIX_def=[LOAD_MATRIX_new(1:2*nx,:);zeros(2*nx*nz-2*nx,3);LOAD_MATRIX_new(2*nx+1:end,:)];
Recovery_Matrix_def=[Recovery_Matrix_new(:,1:2*nx),zeros(size(Recovery_Matrix_new,1),2*nx*nz-2*nx),Recovery_Matrix_new(:,2*nx+1:end)];
FAN_Center_Recovery_vector_def=[FAN_Center_Recovery_vector_new(:,1:2*nx),zeros(size(FAN_Center_Recovery_vector_new,1),2*nx*nz-2*nx),FAN_Center_Recovery_vector_new(:,2*nx+1:end)];
[X,Z]=meshgrid(X_dis,Z_dis);
X=X(:);
Z=Z(:);
[Z,Zid]=sort(Z);
X=X(Zid);
Y=zeros(size(X));
DESIGN_ZONE_COORD=[(1:length(X)).',X,Y,Z];
El_id=(1:(nx-1)*(nz-1)).';
first_node_index=El_id+fix((El_id-1)/(nx-1));
second_node_index=El_id+fix((El_id-1)/(nx-1))+nx;
third_node_index=El_id+fix((El_id-1)/(nx-1))+nx+1;
fourth_node_index=El_id+fix((El_id-1)/(nx-1))+1;
DESIGN_ZONE_ELEMENT_NODES=[first_node_index,second_node_index,third_node_index,fourth_node_index];
first_node_DOFs_index=[2*(first_node_index-1)+1,2*(first_node_index-1)+2];
second_node_DOFs_index=[2*(second_node_index-1)+1,2*(second_node_index-1)+2];
third_node_DOFs_index=[2*(third_node_index-1)+1,2*(third_node_index-1)+2];
fourth_node_DOFs_index=[2*(fourth_node_index-1)+1,2*(fourth_node_index-1)+2];
DESIGN_ZONE_ELEMENT_DOFS=[first_node_DOFs_index,second_node_DOFs_index,third_node_DOFs_index,fourth_node_DOFs_index];
DESIGN_ZONE_ELEMENT_X=reshape(DESIGN_ZONE_COORD(DESIGN_ZONE_ELEMENT_NODES,2),size(DESIGN_ZONE_ELEMENT_NODES));
DESIGN_ZONE_ELEMENT_Z=reshape(DESIGN_ZONE_COORD(DESIGN_ZONE_ELEMENT_NODES,4),size(DESIGN_ZONE_ELEMENT_NODES));
x1=DESIGN_ZONE_ELEMENT_X(:,1);
x2=DESIGN_ZONE_ELEMENT_X(:,2);
x3=DESIGN_ZONE_ELEMENT_X(:,3);
x4=DESIGN_ZONE_ELEMENT_X(:,4);
z1=DESIGN_ZONE_ELEMENT_Z(:,1);
z2=DESIGN_ZONE_ELEMENT_Z(:,2);
z3=DESIGN_ZONE_ELEMENT_Z(:,3);
z4=DESIGN_ZONE_ELEMENT_Z(:,4);
A2=(x4-x2).*(z3-z1)-(x3-x1).*(z4-z2);
E0=210000;
Emin=0.1;
nu=0.3;
%matrerial constant for unitary young modul
C1=(1-nu)/(1+nu)/(1-2*nu); %E1/E
C2=C1*nu/(1-nu); %E2/E
C3=1/2/(1+nu); %G/E
%Transofrmations
T1=sparse(1:8,[2,3,4,1,6,7,8,5],ones(1,8),8,8);
T2=sparse(1:8,[7,6,5,8,3,2,1,4],ones(1,8),8,8);
T3=sparse(1:8,[5,6,7,8,1,2,3,4],ones(1,8),8,8);
%% Group A k11
Estar=C1;
G=C3;
s1=2*(z4-z2).^2;
s2=2*(x4-x2).^2;
s3=-s1/2;
s4=-s2/2;
t1=(z2-z3).^2+(z3-z4).^2+(z4-z2).^2;
t2=(x2-x3).^2+(x3-x4).^2+(x4-x2).^2;
t3=(z4-z3).^2+(z3-z2).^2;
t4=(x4-x3).^2+(x3-x2).^2;
%transformation dipendent terms
f1=(x1+x3).*(z4-z2)-(z1+z3).*(x4-x2)-2*(x2.*z4-x4.*z2);
f2=(z2+z4).*(x3-x1)-(x2+x4).*(z3-z1)-2*(x3.*z1-x1.*z3);
I=DESIGN_ZONE_ELEMENT_DOFS(:,1);
J=I;
k_ij=1/2*((A2.*(Estar*s1+G*s2)+f1.*(Estar*s3+G*s4))./(3*A2.^2-f1.^2)+(A2.*(Estar*t1+G*t2)+f2.*(Estar*t3+G*t4))./(3*A2.^2-f2.^2));
DOF_id=[3,5,7,8,2,4,6];
for dof=DOF_id
    if dof==2
        T=T2;
    else
        T=T1;
    end
    coord_new=(T*[x1,x2,x3,x4,z1,z2,z3,z4].').';
    x1=coord_new(:,1);
    x2=coord_new(:,2);
    x3=coord_new(:,3);
    x4=coord_new(:,4);
    z1=coord_new(:,5);
    z2=coord_new(:,6);
    z3=coord_new(:,7);
    z4=coord_new(:,8);
    %transformation dipendent terms
    f1=(x1+x3).*(z4-z2)-(z1+z3).*(x4-x2)-2*(x2.*z4-x4.*z2);
    f2=(z2+z4).*(x3-x1)-(x2+x4).*(z3-z1)-2*(x3.*z1-x1.*z3);
    Estar=C1;
    G=C3;
    s1=2*(z4-z2).^2;
    s2=2*(x4-x2).^2;
    s3=-s1/2;
    s4=-s2/2;
    t1=(z2-z3).^2+(z3-z4).^2+(z4-z2).^2;
    t2=(x2-x3).^2+(x3-x4).^2+(x4-x2).^2;
    t3=(z4-z3).^2+(z3-z2).^2;
    t4=(x4-x3).^2+(x3-x2).^2;
    I=[I;DESIGN_ZONE_ELEMENT_DOFS(:,dof)];
    J=[J;DESIGN_ZONE_ELEMENT_DOFS(:,dof)];
    k_ij=[k_ij;1/2*((A2.*(Estar*s1+G*s2)+f1.*(Estar*s3+G*s4))./(3*A2.^2-f1.^2)+(A2.*(Estar*t1+G*t2)+f2.*(Estar*t3+G*t4))./(3*A2.^2-f2.^2))];
end
%%
%Group B k21
x1=DESIGN_ZONE_ELEMENT_X(:,1);
x2=DESIGN_ZONE_ELEMENT_X(:,2);
x3=DESIGN_ZONE_ELEMENT_X(:,3);
x4=DESIGN_ZONE_ELEMENT_X(:,4);
z1=DESIGN_ZONE_ELEMENT_Z(:,1);
z2=DESIGN_ZONE_ELEMENT_Z(:,2);
z3=DESIGN_ZONE_ELEMENT_Z(:,3);
z4=DESIGN_ZONE_ELEMENT_Z(:,4);
Estar=C2;
G=C3;
s1=2*(x2-x4).*(z4-z2);
s2=s1;
s3=-s1/2;
s4=3;
t1=x2.*(z4-2*z2+z3)+x3.*(z2-2*z3+z4)+x4.*(z2-2*z4+z3);
t2=t1;
t3=x2.*(z2-z3)+x3.*(z4-z2)+x4.*(z3-z4);
t4=t3;
%transformation dipendent terms
f1=(x1+x3).*(z4-z2)-(z1+z3).*(x4-x2)-2*(x2.*z4-x4.*z2);
f2=(z2+z4).*(x3-x1)-(x2+x4).*(z3-z1)-2*(x3.*z1-x1.*z3);
I=[I;DESIGN_ZONE_ELEMENT_DOFS(:,2)];
J=[J;DESIGN_ZONE_ELEMENT_DOFS(:,1)];
k_ij=[k_ij;1/2*((A2.*(Estar*s1+G*s2)+f1.*(Estar*s3+G*s4))./(3*A2.^2-f1.^2)+(A2.*(Estar*t1+G*t2)+f2.*(Estar*t3+G*t4))./(3*A2.^2-f2.^2))];
DOF_id_i=[4,6,8];
DOF_id_j=[3,5,7];
for df=1:length(DOF_id_i)
    dofi=DOF_id_i(df);
    dofj=DOF_id_j(df);
    T=T1;
    coord_new=(T*[x1,x2,x3,x4,z1,z2,z3,z4].').';
    x1=coord_new(:,1);
    x2=coord_new(:,2);
    x3=coord_new(:,3);
    x4=coord_new(:,4);
    z1=coord_new(:,5);
    z2=coord_new(:,6);
    z3=coord_new(:,7);
    z4=coord_new(:,8);
    %transformation dipendent terms
    f1=(x1+x3).*(z4-z2)-(z1+z3).*(x4-x2)-2*(x2.*z4-x4.*z2);
    f2=(z2+z4).*(x3-x1)-(x2+x4).*(z3-z1)-2*(x3.*z1-x1.*z3);
    Estar=C2;
    G=C3;
    s1=2*(x2-x4).*(z4-z2);
    s2=s1;
    s3=-s1/2;
    s4=3;
    t1=x2.*(z4-2*z2+z3)+x3.*(z2-2*z3+z4)+x4.*(z2-2*z4+z3);
    t2=t1;
    t3=x2.*(z2-z3)+x3.*(z4-z2)+x4.*(z3-z4);
    t4=t3;
    I=[I;DESIGN_ZONE_ELEMENT_DOFS(:,dofi)];
    J=[J;DESIGN_ZONE_ELEMENT_DOFS(:,dofj)];
    k_ij=[k_ij;1/2*((A2.*(Estar*s1+G*s2)+f1.*(Estar*s3+G*s4))./(3*A2.^2-f1.^2)+(A2.*(Estar*t1+G*t2)+f2.*(Estar*t3+G*t4))./(3*A2.^2-f2.^2))];
end
%
%% Group C k31
x1=DESIGN_ZONE_ELEMENT_X(:,1);
x2=DESIGN_ZONE_ELEMENT_X(:,2);
x3=DESIGN_ZONE_ELEMENT_X(:,3);
x4=DESIGN_ZONE_ELEMENT_X(:,4);
z1=DESIGN_ZONE_ELEMENT_Z(:,1);
z2=DESIGN_ZONE_ELEMENT_Z(:,2);
z3=DESIGN_ZONE_ELEMENT_Z(:,3);
z4=DESIGN_ZONE_ELEMENT_Z(:,4);
Estar=C1;
G=C3;
s1=(z4-z2).*(2*z1-z3-z4);
s2=(x4-x2).*(2*x1-x3-x4);
s3=(z4-z2).*(z4-z1);
s4=(x4-x2).*(x4-x1);
t1=(z3-z1).*(2*z2-z3-z4);
t2=(x3-x1).*(2*x2-x3-x4);
t3=(z3-z1).*(z3-z2);
t4=(x3-x1).*(x3-x2);
%transformation dipendent terms
f1=(x1+x3).*(z4-z2)-(z1+z3).*(x4-x2)-2*(x2.*z4-x4.*z2);
f2=(z2+z4).*(x3-x1)-(x2+x4).*(z3-z1)-2*(x3.*z1-x1.*z3);
I=[I;DESIGN_ZONE_ELEMENT_DOFS(:,3)];
J=[J;DESIGN_ZONE_ELEMENT_DOFS(:,1)];
k_ij=[k_ij;1/2*((A2.*(Estar*s1+G*s2)+f1.*(Estar*s3+G*s4))./(3*A2.^2-f1.^2)+(A2.*(Estar*t1+G*t2)+f2.*(Estar*t3+G*t4))./(3*A2.^2-f2.^2))];
DOF_id_i=[5,7,7,8,8,4,6];
DOF_id_j=[3,5,1,6,2,2,4];
for df=1:length(DOF_id_i)
    dofi=DOF_id_i(df);
    dofj=DOF_id_j(df);
    if dofi==8&&dofj==6
        T=T2;
    else
        T=T1;
    end
    coord_new=(T*[x1,x2,x3,x4,z1,z2,z3,z4].').';
    x1=coord_new(:,1);
    x2=coord_new(:,2);
    x3=coord_new(:,3);
    x4=coord_new(:,4);
    z1=coord_new(:,5);
    z2=coord_new(:,6);
    z3=coord_new(:,7);
    z4=coord_new(:,8);
    %transformation dipendent terms
    f1=(x1+x3).*(z4-z2)-(z1+z3).*(x4-x2)-2*(x2.*z4-x4.*z2);
    f2=(z2+z4).*(x3-x1)-(x2+x4).*(z3-z1)-2*(x3.*z1-x1.*z3);
    Estar=C1;
    G=C3;
    s1=(z4-z2).*(2*z1-z3-z4);
    s2=(x4-x2).*(2*x1-x3-x4);
    s3=(z4-z2).*(z4-z1);
    s4=(x4-x2).*(x4-x1);
    t1=(z3-z1).*(2*z2-z3-z4);
    t2=(x3-x1).*(2*x2-x3-x4);
    t3=(z3-z1).*(z3-z2);
    t4=(x3-x1).*(x3-x2);
    I=[I;DESIGN_ZONE_ELEMENT_DOFS(:,dofi)];
    J=[J;DESIGN_ZONE_ELEMENT_DOFS(:,dofj)];
    k_ij=[k_ij;1/2*((A2.*(Estar*s1+G*s2)+f1.*(Estar*s3+G*s4))./(3*A2.^2-f1.^2)+(A2.*(Estar*t1+G*t2)+f2.*(Estar*t3+G*t4))./(3*A2.^2-f2.^2))];
end
%% Group D k41
x1=DESIGN_ZONE_ELEMENT_X(:,1);
x2=DESIGN_ZONE_ELEMENT_X(:,2);
x3=DESIGN_ZONE_ELEMENT_X(:,3);
x4=DESIGN_ZONE_ELEMENT_X(:,4);
z1=DESIGN_ZONE_ELEMENT_Z(:,1);
z2=DESIGN_ZONE_ELEMENT_Z(:,2);
z3=DESIGN_ZONE_ELEMENT_Z(:,3);
z4=DESIGN_ZONE_ELEMENT_Z(:,4);
Estar=C2;
G=C3;
s1=(x3-x1).*(z4-z2)+(x4-x1).*(z4-z2);
s2=(z3-z1).*(x4-x2)+(z4-z1).*(x4-x2);
s3=(x4-x1).*(z2-z4);
s4=(z4-z1).*(x2-x4);
t1=(x3-x1).*(z4-z2)+(x3-x1).*(z3-z2);
t2=(z3-z1).*(x4-x2)+(z3-z1).*(x3-x2);
t3=(x3-x1).*(z2-z3);
t4=(z3-z1).*(x2-x3);
%transformation dipendent terms
f1=(x1+x3).*(z4-z2)-(z1+z3).*(x4-x2)-2*(x2.*z4-x4.*z2);
f2=(z2+z4).*(x3-x1)-(x2+x4).*(z3-z1)-2*(x3.*z1-x1.*z3);
I=[I;DESIGN_ZONE_ELEMENT_DOFS(:,4)];
J=[J;DESIGN_ZONE_ELEMENT_DOFS(:,1)];
k_ij=[k_ij;1/2*((A2.*(Estar*s1+G*s2)+f1.*(Estar*s3+G*s4))./(3*A2.^2-f1.^2)+(A2.*(Estar*t1+G*t2)+f2.*(Estar*t3+G*t4))./(3*A2.^2-f2.^2))];
DOF_id_i=[6,8,7,3,5,7,8];
DOF_id_j=[3,5,2,2,4,6,1];
for df=1:length(DOF_id_i)
    dofi=DOF_id_i(df);
    dofj=DOF_id_j(df);
    if dofi==3&&dofj==2
        T=T3;
    else
        T=T1;
    end
    coord_new=(T*[x1,x2,x3,x4,z1,z2,z3,z4].').';
    x1=coord_new(:,1);
    x2=coord_new(:,2);
    x3=coord_new(:,3);
    x4=coord_new(:,4);
    z1=coord_new(:,5);
    z2=coord_new(:,6);
    z3=coord_new(:,7);
    z4=coord_new(:,8);
    %transformation dipendent terms
    f1=(x1+x3).*(z4-z2)-(z1+z3).*(x4-x2)-2*(x2.*z4-x4.*z2);
    f2=(z2+z4).*(x3-x1)-(x2+x4).*(z3-z1)-2*(x3.*z1-x1.*z3);
    Estar=C2;
    G=C3;
    s1=(x3-x1).*(z4-z2)+(x4-x1).*(z4-z2);
    s2=(z3-z1).*(x4-x2)+(z4-z1).*(x4-x2);
    s3=(x4-x1).*(z2-z4);
    s4=(z4-z1).*(x2-x4);
    t1=(x3-x1).*(z4-z2)+(x3-x1).*(z3-z2);
    t2=(z3-z1).*(x4-x2)+(z3-z1).*(x3-x2);
    t3=(x3-x1).*(z2-z3);
    t4=(z3-z1).*(x2-x3);
    I=[I;DESIGN_ZONE_ELEMENT_DOFS(:,dofi)];
    J=[J;DESIGN_ZONE_ELEMENT_DOFS(:,dofj)];
    k_ij=[k_ij;1/2*((A2.*(Estar*s1+G*s2)+f1.*(Estar*s3+G*s4))./(3*A2.^2-f1.^2)+(A2.*(Estar*t1+G*t2)+f2.*(Estar*t3+G*t4))./(3*A2.^2-f2.^2))];
end
%% Group E k51
x1=DESIGN_ZONE_ELEMENT_X(:,1);
x2=DESIGN_ZONE_ELEMENT_X(:,2);
x3=DESIGN_ZONE_ELEMENT_X(:,3);
x4=DESIGN_ZONE_ELEMENT_X(:,4);
z1=DESIGN_ZONE_ELEMENT_Z(:,1);
z2=DESIGN_ZONE_ELEMENT_Z(:,2);
z3=DESIGN_ZONE_ELEMENT_Z(:,3);
z4=DESIGN_ZONE_ELEMENT_Z(:,4);
Estar=C1;
G=C3;
s1=-(z4-z2).^2;
s2=-(x4-x2).^2;
s3=zeros(size(s1));
s4=zeros(size(s1));
t1=(z3+z1).*(z4+z2)-2*(z4-z2).^2-2*(z1.*z3+z2.*z4);
t2=(x3+x1).*(x4+x2)-2*(x4-x2).^2-2*(x1.*x3+x2.*x4);
t3=(z4-z2).*(z1-z2+z3-z4);
t4=(x4-x2).*(x1-x2+x3-x4);
%transformation dipendent terms
f1=(x1+x3).*(z4-z2)-(z1+z3).*(x4-x2)-2*(x2.*z4-x4.*z2);
f2=(z2+z4).*(x3-x1)-(x2+x4).*(z3-z1)-2*(x3.*z1-x1.*z3);
I=[I;DESIGN_ZONE_ELEMENT_DOFS(:,5)];
J=[J;DESIGN_ZONE_ELEMENT_DOFS(:,1)];
k_ij=[k_ij;1/2*((A2.*(Estar*s1+G*s2)+f1.*(Estar*s3+G*s4))./(3*A2.^2-f1.^2)+(A2.*(Estar*t1+G*t2)+f2.*(Estar*t3+G*t4))./(3*A2.^2-f2.^2))];
DOF_id_i=[7,8,6];
DOF_id_j=[3,4,2];
for df=1:length(DOF_id_i)
    dofi=DOF_id_i(df);
    dofj=DOF_id_j(df);
    if dofi==8&&dofj==4
        T=T2;
    else
        T=T1;
    end
    coord_new=(T*[x1,x2,x3,x4,z1,z2,z3,z4].').';
    x1=coord_new(:,1);
    x2=coord_new(:,2);
    x3=coord_new(:,3);
    x4=coord_new(:,4);
    z1=coord_new(:,5);
    z2=coord_new(:,6);
    z3=coord_new(:,7);
    z4=coord_new(:,8);
    %transformation dipendent terms
    f1=(x1+x3).*(z4-z2)-(z1+z3).*(x4-x2)-2*(x2.*z4-x4.*z2);
    f2=(z2+z4).*(x3-x1)-(x2+x4).*(z3-z1)-2*(x3.*z1-x1.*z3);
    Estar=C1;
    G=C3;
    s1=-(z4-z2).^2;
    s2=-(x4-x2).^2;
    s3=zeros(size(s1));
    s4=zeros(size(s1));
    t1=(z3+z1).*(z4+z2)-2*(z4-z2).^2-2*(z1.*z3+z2.*z4);
    t2=(x3+x1).*(x4+x2)-2*(x4-x2).^2-2*(x1.*x3+x2.*x4);
    t3=(z4-z2).*(z1-z2+z3-z4);
    t4=(x4-x2).*(x1-x2+x3-x4);
    I=[I;DESIGN_ZONE_ELEMENT_DOFS(:,dofi)];
    J=[J;DESIGN_ZONE_ELEMENT_DOFS(:,dofj)];
    k_ij=[k_ij;1/2*((A2.*(Estar*s1+G*s2)+f1.*(Estar*s3+G*s4))./(3*A2.^2-f1.^2)+(A2.*(Estar*t1+G*t2)+f2.*(Estar*t3+G*t4))./(3*A2.^2-f2.^2))];
end
%% Group F k61
x1=DESIGN_ZONE_ELEMENT_X(:,1);
x2=DESIGN_ZONE_ELEMENT_X(:,2);
x3=DESIGN_ZONE_ELEMENT_X(:,3);
x4=DESIGN_ZONE_ELEMENT_X(:,4);
z1=DESIGN_ZONE_ELEMENT_Z(:,1);
z2=DESIGN_ZONE_ELEMENT_Z(:,2);
z3=DESIGN_ZONE_ELEMENT_Z(:,3);
z4=DESIGN_ZONE_ELEMENT_Z(:,4);
Estar=C2;
G=C3;
s1=(x4-x2).*(z4-z2);
s2=s1;
s3=zeros(size(s1));
s4=zeros(size(s1));
t1=(x4-x2).*(z4-z2)+(x2-x1).*(z2-z3)+(x4-x1).*(z4-z3);
t2=(z4-z2).*(x4-x2)+(z2-z1).*(x2-x3)+(z4-z1).*(x4-x3);
t3=(x2-x1).*(z3-z2)+(x4-x1).*(z4-z3);
t4=(x2-x1).*(z3-z2)+(x4-x1).*(z4-z3);
%transformation dipendent terms
f1=(x1+x3).*(z4-z2)-(z1+z3).*(x4-x2)-2*(x2.*z4-x4.*z2);
f2=(z2+z4).*(x3-x1)-(x2+x4).*(z3-z1)-2*(x3.*z1-x1.*z3);
I=[I;DESIGN_ZONE_ELEMENT_DOFS(:,6)];
J=[J;DESIGN_ZONE_ELEMENT_DOFS(:,1)];
k_ij=[k_ij;1/2*((A2.*(Estar*s1+G*s2)+f1.*(Estar*s3+G*s4))./(3*A2.^2-f1.^2)+(A2.*(Estar*t1+G*t2)+f2.*(Estar*t3+G*t4))./(3*A2.^2-f2.^2))];
DOF_id_i=[8,5,7];
DOF_id_j=[3,2,4];
for df=1:length(DOF_id_i)
    dofi=DOF_id_i(df);
    dofj=DOF_id_j(df);
    if dofi==5&&dofj==2
        T=T3;
    else
        T=T1;
    end
    coord_new=(T*[x1,x2,x3,x4,z1,z2,z3,z4].').';
    x1=coord_new(:,1);
    x2=coord_new(:,2);
    x3=coord_new(:,3);
    x4=coord_new(:,4);
    z1=coord_new(:,5);
    z2=coord_new(:,6);
    z3=coord_new(:,7);
    z4=coord_new(:,8);
    %transformation dipendent terms
    f1=(x1+x3).*(z4-z2)-(z1+z3).*(x4-x2)-2*(x2.*z4-x4.*z2);
    f2=(z2+z4).*(x3-x1)-(x2+x4).*(z3-z1)-2*(x3.*z1-x1.*z3);
    Estar=C2;
    G=C3;
    s1=(x4-x2).*(z4-z2);
    s2=s1;
    s3=zeros(size(s1));
    s4=zeros(size(s1));
    t1=(x4-x2).*(z4-z2)+(x2-x1).*(z2-z3)+(x4-x1).*(z4-z3);
    t2=(z4-z2).*(x4-x2)+(z2-z1).*(x2-x3)+(z4-z1).*(x4-x3);
    t3=(x2-x1).*(z3-z2)+(x4-x1).*(z4-z3);
    t4=(x2-x1).*(z3-z2)+(x4-x1).*(z4-z3);
    I=[I;DESIGN_ZONE_ELEMENT_DOFS(:,dofi)];
    J=[J;DESIGN_ZONE_ELEMENT_DOFS(:,dofj)];
    k_ij=[k_ij;1/2*((A2.*(Estar*s1+G*s2)+f1.*(Estar*s3+G*s4))./(3*A2.^2-f1.^2)+(A2.*(Estar*t1+G*t2)+f2.*(Estar*t3+G*t4))./(3*A2.^2-f2.^2))];
end
% BCs
% fixed DOFs in the region x>=4000&Z==max(Z)
fixed_nodes=DESIGN_ZONE_COORD(DESIGN_ZONE_COORD(:,2)>=4000&DESIGN_ZONE_COORD(:,4)==max(DESIGN_ZONE_COORD(:,4)),1);
fixed_dofs=[2*(fixed_nodes-1)+1;2*(fixed_nodes-1)+2];
fixed_dofs=sort(fixed_dofs);
All_DOFs=(1:2*nx*nz).';
free_dofs=setdiff(All_DOFs,fixed_dofs);
Lagrangian_multip=new_DOFs_MAP(end-sum(new_DOFs_MAP(:,1)==0)+1:end,4);
free_dofs=union(free_dofs,Lagrangian_multip);
All_DOFs=union(All_DOFs,Lagrangian_multip);
U=zeros(length(All_DOFs),1);

%% Preparing Filter
rmin=1;
nelx=nx-1;
nely=nz-1;
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
volfrac=0.4; penal=3; ft=1;
x = repmat(volfrac,nely,nelx);
xPhys = x;
loop = 0;
change = 1;
m = 2;
n = length(xPhys(:));
epsimin = 0.0000001;
eeen    = ones(n,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = xPhys(:);
xold1   = xval;
xold2   = xval;
xmin    = zeron;
xmax    = eeen;
low     = xmin;
upp     = xmax;
c       = 1000*eeem;
d       = eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 100;
kkttol  = 0;
load('gamma_structure')
load('U_tilder')
%
%%%% The iterations start:
kktnorm = kkttol+10;
outit = 0;
%% FE-ANALYSIS
E=Emin+xPhys(:).^penal*(E0-Emin);
E=repmat(E,length(k_ij)/length(E),1);
dK=penal*(E0-Emin)*xPhys(:).^(penal-1);
dK=repmat(dK,length(k_ij)/length(dK),1);
Lind=(1:length(xPhys(:))).';
Lind=repmat(Lind,length(k_ij)/length(Lind),1);
dK=k_ij.*dK;
sK=k_ij.*E;
K = sparse(I,J,sK,2*nx*nz+sum(DOFs_MAP(:,1)==0),2*nx*nz+sum(DOFs_MAP(:,1)==0)); K=K+triu(K,1);
K=K_def+K;
F=[LOAD_MATRIX_def(:,1),Recovery_Matrix_def.',FAN_Center_Recovery_vector_def.'];
Res=K(free_dofs,free_dofs)\F(free_dofs,:);
U(free_dofs)=Res(:,1);
dK_dxU=zeros(length(U),length(xPhys(:)));
for xdim=1:length(xPhys(:))
    dK_dx=sparse(I,J,dK.*(Lind==xdim),2*nx*nz+sum(DOFs_MAP(:,1)==0),2*nx*nz+sum(DOFs_MAP(:,1)==0)); dK_dx=dK_dx+triu(dK_dx,1);
    dK_dxU(:,xdim)=dK_dx*U;
end
H_fan_disp=zeros(size(U,1),1);
H_fan_disp(free_dofs,:)=Res(:,end);
H_s=zeros(size(U,1),size(Res,2)-2);
H_s(free_dofs,:)=Res(:,2:end-1);
%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
Utip=Recovery_Matrix_def*U;
perfo=0;
psi_tild=0;
for s=1:20
    tip_vector=Gamma(s).mat*(U_static(:,1)+Utip);
    R=rms(tip_vector);
    perfo=perfo+Gamma(s).lambda*R;
    N=length(tip_vector);
    R=ones(N,1)*R;
    psi_tild=psi_tild-Gamma(s).lambda/N*H_s*Gamma(s).mat.'*(tip_vector./R);
end
FAN_Ax_disp=U_FAN_static(1)+FAN_Center_Recovery_vector_def*U;
dperfo=psi_tild.'*dK_dxU;
dfandisp=-H_fan_disp.'*dK_dxU;
dv = A2/2;
%% FILTERING/MODIFICATION OF SENSITIVITIES
if ft == 1
    dperfo(:) = H*(x(:).*dperfo.')./Hs./max(1e-3,x(:));
elseif ft == 2
    dperfo = H*(dperfo.'./Hs);
    dv(:) = H*(dv(:)./Hs);
end
f0val=perfo;
fval=[sum(xPhys(:).*A2/2) - volfrac*sum(A2/2);(-0.01-FAN_Ax_disp)/0.01];
df0dx=reshape(dperfo,[],1);
dfdx=[dv(:)';-dfandisp/0.01];
innerit=0;
outvector1 = [outeriter innerit f0val fval'];
outvector2 = xval;
h = figure(1); set(h,'Color',[1 1 1]);
hold on; h24_patch = patch('Vertices',DESIGN_ZONE_COORD(:,[2,4]),'Faces',DESIGN_ZONE_ELEMENT_NODES,'FaceVertexCData',(1-xval)*[1 1 1],'FaceColor','flat'); axis equal; axis off; drawnow;hold off
fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',outeriter,perfo, ...
    mean(xPhys(:)),kktnorm);
figure(2)
hold on
scatter(outeriter,f0val,'fill','k')
figure(3)
hold on
scatter(outeriter,FAN_Ax_disp,'fill','b')

%% START ITERATION
while kktnorm > kkttol & outit < maxoutit
    outit   = outit+1;
    outeriter = outeriter+1;
    %% MMA code optimization
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;
    if ft == 1
        xPhys=xval;
    elseif ft == 2
        xPhys = (H*xval(:))./Hs;
    end
    %% FE-ANALYSIS
    E=Emin+xPhys(:).^penal*(E0-Emin);
    E=repmat(E,length(k_ij)/length(E),1);
    dK=penal*(E0-Emin)*xPhys(:).^(penal-1);
    dK=repmat(dK,length(k_ij)/length(dK),1);
    Lind=(1:length(xPhys(:))).';
    Lind=repmat(Lind,length(k_ij)/length(Lind),1);
    dK=k_ij.*dK;
    sK=k_ij.*E;
    K = sparse(I,J,sK,2*nx*nz+sum(DOFs_MAP(:,1)==0),2*nx*nz+sum(DOFs_MAP(:,1)==0)); K=K+triu(K,1);
    K=K_def+K;
    F=[LOAD_MATRIX_def(:,1),Recovery_Matrix_def.',FAN_Center_Recovery_vector_def.'];
    Res=K(free_dofs,free_dofs)\F(free_dofs,:);
    U(free_dofs)=Res(:,1);
    dK_dxU=zeros(length(U),length(xPhys(:)));
    for xdim=1:length(xPhys(:))
        dK_dx=sparse(I,J,dK.*(Lind==xdim),2*nx*nz+sum(DOFs_MAP(:,1)==0),2*nx*nz+sum(DOFs_MAP(:,1)==0)); dK_dx=dK_dx+triu(dK_dx,1);
        dK_dxU(:,xdim)=dK_dx*U;
    end
    H_fan_disp=zeros(size(U,1),1);
    H_fan_disp(free_dofs,:)=Res(:,end);
    H_s=zeros(size(U,1),size(Res,2)-2);
    H_s(free_dofs,:)=Res(:,2:end-1);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    Utip=Recovery_Matrix_def*U;
    perfo=0;
    psi_tild=0;
    for s=1:20
        tip_vector=Gamma(s).mat*(U_static(:,1)+Utip);
        R=rms(tip_vector);
        perfo=perfo+Gamma(s).lambda*R;
        N=length(tip_vector);
        R=ones(N,1)*R;
        psi_tild=psi_tild-Gamma(s).lambda/N*H_s*Gamma(s).mat.'*(tip_vector./R);
    end
    FAN_Ax_disp=U_FAN_static(1)+FAN_Center_Recovery_vector_def*U;
    dperfo=psi_tild.'*dK_dxU;
    dfandisp=-H_fan_disp.'*dK_dxU;
    dv = A2/2;
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dperfo(:) = H*(x(:).*dperfo.')./Hs./max(1e-3,x(:));
    elseif ft == 2
        dperfo = H*(dperfo.'./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
    f0val=perfo;
    fval=[sum(xPhys(:).*A2/2) - volfrac*sum(A2/2);(-0.01-FAN_Ax_disp)/0.01];
    df0dx=reshape(dperfo,[],1);
    dfdx=[dv(:)';-dfandisp/0.01];
    
    %% %% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
    outvector1 = [outeriter innerit f0val fval'];
    outvector2 = xval';
    xPhys=xval';
    %% PLOT DENSITIES
    
    %     colormap(gray);h = figure(1); set(h,'Color',[1 1 1]);
    h = figure(1); set(h,'Color',[1 1 1]);
    hold on; h24_patch = patch('Vertices',DESIGN_ZONE_COORD(:,[2,4]),'Faces',DESIGN_ZONE_ELEMENT_NODES,'FaceVertexCData',(1-xval)*[1 1 1],'FaceColor','flat'); axis equal; axis off; drawnow;hold off
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',outeriter,perfo, ...
        sum(xPhys(:).*A2/2)/sum(A2/2),kktnorm);
    figure(2)
    hold on
    scatter(outeriter,f0val,'fill','k')
    figure(3)
    hold on
    scatter(outeriter,FAN_Ax_disp,'fill','b')
end