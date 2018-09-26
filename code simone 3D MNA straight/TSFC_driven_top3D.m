clear
close all
load('substructure_results')
load('st3D.dat.mat')
load('dz3Dstr.mat')
load('connectivity_interface')
load('gamma_structure')
load('U_tilder')
%the point at the same x have the same 1 and 3 diplacement + linear
%% interpolation.
%find the interpolation nodes of the Design ZONE
DESIGN_ZONE_COORD=[(1:size(COORD,1)).',COORD];
El_id=(1:size(ELEMENT,1)).';
first_node_index=ELEMENT(:,2);
second_node_index=ELEMENT(:,3);
third_node_index=ELEMENT(:,4);
fourth_node_index=ELEMENT(:,5);
fivth_node_index=ELEMENT(:,6);
sixth_node_index=ELEMENT(:,7);
seventh_node_index=ELEMENT(:,8);
eigth_node_index=ELEMENT(:,9);
DESIGN_ZONE_ELEMENT_NODES=[first_node_index,second_node_index,third_node_index,fourth_node_index,fivth_node_index,sixth_node_index,seventh_node_index,eigth_node_index];
first_node_DOFs_index=[3*(first_node_index-1)+1,3*(first_node_index-1)+2,3*(first_node_index-1)+3];
second_node_DOFs_index=[3*(second_node_index-1)+1,3*(second_node_index-1)+2,3*(second_node_index-1)+3];
third_node_DOFs_index=[3*(third_node_index-1)+1,3*(third_node_index-1)+2,3*(third_node_index-1)+3];
fourth_node_DOFs_index=[3*(fourth_node_index-1)+1,3*(fourth_node_index-1)+2,3*(fourth_node_index-1)+3];
fivth_node_DOFs_index=[3*(fivth_node_index-1)+1,3*(fivth_node_index-1)+2,3*(fivth_node_index-1)+3];
sixth_node_DOFs_index=[3*(sixth_node_index-1)+1,3*(sixth_node_index-1)+2,3*(sixth_node_index-1)+3];
seventh_node_DOFs_index=[3*(seventh_node_index-1)+1,3*(seventh_node_index-1)+2,3*(seventh_node_index-1)+3];
eighth_node_DOFs_index=[3*(eigth_node_index-1)+1,3*(eigth_node_index-1)+2,3*(eigth_node_index-1)+3];
DESIGN_ZONE_ELEMENT_DOFS=[first_node_DOFs_index,second_node_DOFs_index,third_node_DOFs_index,fourth_node_DOFs_index,fivth_node_DOFs_index,sixth_node_DOFs_index,seventh_node_DOFs_index,eighth_node_DOFs_index];
new_DOFs_MAP=zeros(3*size(COORD,1)+sum(DOFs_MAP(:,1)==0),4);
for k=1:size(COORD,1)
    new_DOFs_MAP(3*k-2,:)=[k,k,1,3*k-2];
    new_DOFs_MAP(3*k-1,:)=[k,k,2,3*k-1];
    new_DOFs_MAP(3*k,:)=[k,k,3,3*k];
end
new_DOFs_MAP(end-sum(DOFs_MAP(:,1)==0)+1:end,:)=[DOFs_MAP(end-sum(DOFs_MAP(:,1)==0)+1:end,1:3),(3*k+(1:sum(DOFs_MAP(:,1)==0)))'];
interface_DZ=NODE_SET.NID;
interface_DZ_coord=COORD(interface_DZ,:);
[ID4]=ismember(FACES(:,[2:5]),interface_DZ);
IDF=ID4(:,1)&ID4(:,2)&ID4(:,3)&ID4(:,4);
FACE_intdz=FACES(IDF,:);
x1dz=COORD(FACE_intdz(:,2),1);
y1dz=COORD(FACE_intdz(:,2),2);
z1dz=COORD(FACE_intdz(:,2),3);
x2dz=COORD(FACE_intdz(:,3),1);
y2dz=COORD(FACE_intdz(:,3),2);
z2dz=COORD(FACE_intdz(:,3),3);
x3dz=COORD(FACE_intdz(:,4),1);
y3dz=COORD(FACE_intdz(:,4),2);
z3dz=COORD(FACE_intdz(:,4),3);
x4dz=COORD(FACE_intdz(:,5),1);
y4dz=COORD(FACE_intdz(:,5),2);
z4dz=COORD(FACE_intdz(:,5),3);
%interface Mass Matrix

N1=@(csi,eta) 1/4*(1-csi).*(1-eta);
N2=@(csi,eta) 1/4*(1-csi).*(1+eta);
N3=@(csi,eta) 1/4*(1+csi).*(1+eta);
N4=@(csi,eta) 1/4*(1+csi).*(1-eta);
%Interface Mass Matrix calculation
%Lagrange function on the interface:
dN1_dcsi=@(csi,eta) -1/4*(1-eta);
dN1_deta=@(csi,eta) -1/4*(1-csi);
dN2_dcsi=@(csi,eta) -1/4*(1+eta);
dN2_deta=@(csi,eta) +1/4*(1-csi);
dN3_dcsi=@(csi,eta) +1/4*(1+eta);
dN3_deta=@(csi,eta) +1/4*(1+csi);
dN4_dcsi=@(csi,eta) +1/4*(1-eta);
dN4_deta=@(csi,eta) -1/4*(1+csi);
Csi=[-1/sqrt(3),-1/sqrt(3),1/sqrt(3),1/sqrt(3)];
Eta=[-1/sqrt(3),1/sqrt(3),1/sqrt(3),-1/sqrt(3)];
for l=1:4
    cs=Csi(l);
    et=Eta(l);
    N1pg(:,l)=N1(cs,et)*ones(size(x1dz));
    N2pg(:,l)=N2(cs,et)*ones(size(x1dz));
    N3pg(:,l)=N3(cs,et)*ones(size(x1dz));
    N4pg(:,l)=N4(cs,et)*ones(size(x1dz));
    dN1_dCsi(:,l)=dN1_dcsi(cs,et);
    dN2_dCsi(:,l)=dN2_dcsi(cs,et);
    dN3_dCsi(:,l)=dN3_dcsi(cs,et);
    dN4_dCsi(:,l)=dN4_dcsi(cs,et);
    dN1_dEta(:,l)=dN1_deta(cs,et);
    dN2_dEta(:,l)=dN2_deta(cs,et);
    dN3_dEta(:,l)=dN3_deta(cs,et);
    dN4_dEta(:,l)=dN4_deta(cs,et);
    dx_dcsi(:,l)=x1dz*dN1_dcsi(cs,et)+x2dz*dN2_dcsi(cs,et)+x3dz*dN3_dcsi(cs,et)+x4dz*dN4_dcsi(cs,et);
    dx_deta(:,l)=x1dz*dN1_deta(cs,et)+x2dz*dN2_deta(cs,et)+x3dz*dN3_deta(cs,et)+x4dz*dN4_deta(cs,et);
    dy_dcsi(:,l)=y1dz*dN1_dcsi(cs,et)+y2dz*dN2_dcsi(cs,et)+y3dz*dN3_dcsi(cs,et)+y4dz*dN4_dcsi(cs,et);
    dy_deta(:,l)=y1dz*dN1_deta(cs,et)+y2dz*dN2_deta(cs,et)+y3dz*dN3_deta(cs,et)+y4dz*dN4_deta(cs,et);
    dz_dcsi(:,l)=z1dz*dN1_dcsi(cs,et)+z2dz*dN2_dcsi(cs,et)+z3dz*dN3_dcsi(cs,et)+z4dz*dN4_dcsi(cs,et);
    dz_deta(:,l)=z1dz*dN1_deta(cs,et)+z2dz*dN2_deta(cs,et)+z3dz*dN3_deta(cs,et)+z4dz*dN4_deta(cs,et);
end
detJdz=sqrt((dy_dcsi.*dz_deta-dy_deta.*dz_dcsi).^2+(-dx_dcsi.*dz_deta+dx_deta.*dz_dcsi).^2+(dx_dcsi.*dy_deta-dx_deta.*dy_dcsi).^2);
DESIGN_ZONE_FACES_DOFS=[3*(FACE_intdz(:,2)-1)+1,3*(FACE_intdz(:,2)-1)+2,3*(FACE_intdz(:,2)-1)+3,3*(FACE_intdz(:,3)-1)+1,3*(FACE_intdz(:,3)-1)+2,3*(FACE_intdz(:,3)-1)+3,3*(FACE_intdz(:,4)-1)+1,3*(FACE_intdz(:,4)-1)+2,3*(FACE_intdz(:,4)-1)+3,3*(FACE_intdz(:,5)-1)+1,3*(FACE_intdz(:,5)-1)+2,3*(FACE_intdz(:,5)-1)+3];
M_e=zeros(12,12,size(detJdz,1));
% I=zeros(30*size(detJdz,1),1);
% J=zeros(30*size(detJdz,1),1);
I=zeros(78*size(detJdz,1),1);
J=zeros(78*size(detJdz,1),1);
m_ij=J;
el_ij=m_ij;
BTB=zeros(12,12,size(detJdz,1),4);
LEN_tot=0;
for l=1:size(detJdz,1)
    for m=1:4
        B=[N1pg(l,m)*eye(3),N2pg(l,m)*eye(3),N3pg(l,m)*eye(3),N4pg(l,m)*eye(3)];
        BTB(:,:,l,m)=B.'*B;
    end
    M_e(:,:,l)=detJdz(l,1)*BTB(:,:,l,1)+detJdz(l,2)*BTB(:,:,l,2)+detJdz(l,3)*BTB(:,:,l,3)+detJdz(l,4)*BTB(:,:,l,4);
    [IND_I,IND_J,Mv]=find(tril(M_e(:,:,l)));
    LEN=length(IND_I);
    I(LEN_tot+(1:LEN)',1)=reshape(DESIGN_ZONE_FACES_DOFS(l,IND_I),[],1);
    J(LEN_tot+(1:LEN)',1)=reshape(DESIGN_ZONE_FACES_DOFS(l,IND_J),[],1);
    m_ij(LEN_tot+(1:LEN)',1)=Mv(:);
    el_ij(LEN_tot+(1:LEN)',1)=l;
    LEN_tot=LEN_tot+LEN;
end
I=I(1:LEN_tot);
J=J(1:LEN_tot);
m_ij=m_ij(1:LEN_tot);
el_ij=el_ij(1:LEN_tot);
i1=I;
j1=J;
i1(I<J)=J(I<J);
j1(I<J)=I(I<J);
I=i1;
J=j1;
MDZ=sparse(I,J,m_ij,size(new_DOFs_MAP,1),size(new_DOFs_MAP,1));MDZ=MDZ+triu(MDZ.',1);MDZ=(MDZ+MDZ.')/2;
[In,Jn]=find(diag(diag(MDZ)));
MDZ=MDZ+(speye(size(MDZ))+sparse(In,Jn,-1*ones(size(In)),3*size(COORD,1)+sum(DOFs_MAP(:,1)==0),3*size(COORD,1)+sum(DOFs_MAP(:,1)==0)));
%same for engine
interface_E=intercoord(:,1);
interface_E_coord=intercoord(:,2:4);
FACE_intE=ELint;
for k=1:size(ELint,1)
    for l=1:size(ELint,2)-1
        FACE_intE(k,l+1)=find(intercoord(:,1)==ELint(k,l+1));
    end
end
ENGINE_FACES_DOFS=[3*(FACE_intE(:,2)-1)+1,3*(FACE_intE(:,2)-1)+2,3*(FACE_intE(:,2)-1)+3,3*(FACE_intE(:,3)-1)+1,3*(FACE_intE(:,3)-1)+2,3*(FACE_intE(:,3)-1)+3,3*(FACE_intE(:,4)-1)+1,3*(FACE_intE(:,4)-1)+2,3*(FACE_intE(:,4)-1)+3,3*(FACE_intE(:,5)-1)+1,3*(FACE_intE(:,5)-1)+2,3*(FACE_intE(:,5)-1)+3];
x1E=interface_E_coord(FACE_intE(:,2),1);
y1E=interface_E_coord(FACE_intE(:,2),2);
z1E=interface_E_coord(FACE_intE(:,2),3);
x2E=interface_E_coord(FACE_intE(:,3),1);
y2E=interface_E_coord(FACE_intE(:,3),2);
z2E=interface_E_coord(FACE_intE(:,3),3);
x3E=interface_E_coord(FACE_intE(:,4),1);
y3E=interface_E_coord(FACE_intE(:,4),2);
z3E=interface_E_coord(FACE_intE(:,4),3);
x4E=interface_E_coord(FACE_intE(:,5),1);
y4E=interface_E_coord(FACE_intE(:,5),2);
z4E=interface_E_coord(FACE_intE(:,5),3);
N1pg=zeros(size(x1E,1),4);
N2pg=N1pg;
N3pg=N1pg;
N4pg=N1pg;
dx_dcsi=N1pg;
dx_deta=N1pg;
dy_dcsi=N1pg;
dy_deta=N1pg;
dz_dcsi=N1pg;
dz_deta=N1pg;
for l=1:4
    cs=Csi(l);
    et=Eta(l);
    N1pg(:,l)=N1(cs,et)*ones(size(x1E));
    N2pg(:,l)=N2(cs,et)*ones(size(x1E));
    N3pg(:,l)=N3(cs,et)*ones(size(x1E));
    N4pg(:,l)=N4(cs,et)*ones(size(x1E));
    dN1_dCsi(:,l)=dN1_dcsi(cs,et);
    dN2_dCsi(:,l)=dN2_dcsi(cs,et);
    dN3_dCsi(:,l)=dN3_dcsi(cs,et);
    dN4_dCsi(:,l)=dN4_dcsi(cs,et);
    dN1_dEta(:,l)=dN1_deta(cs,et);
    dN2_dEta(:,l)=dN2_deta(cs,et);
    dN3_dEta(:,l)=dN3_deta(cs,et);
    dN4_dEta(:,l)=dN4_deta(cs,et);
    dx_dcsi(:,l)=x1E*dN1_dcsi(cs,et)+x2E*dN2_dcsi(cs,et)+x3E*dN3_dcsi(cs,et)+x4E*dN4_dcsi(cs,et);
    dx_deta(:,l)=x1E*dN1_deta(cs,et)+x2E*dN2_deta(cs,et)+x3E*dN3_deta(cs,et)+x4E*dN4_deta(cs,et);
    dy_dcsi(:,l)=y1E*dN1_dcsi(cs,et)+y2E*dN2_dcsi(cs,et)+y3E*dN3_dcsi(cs,et)+y4E*dN4_dcsi(cs,et);
    dy_deta(:,l)=y1E*dN1_deta(cs,et)+y2E*dN2_deta(cs,et)+y3E*dN3_deta(cs,et)+y4E*dN4_deta(cs,et);
    dz_dcsi(:,l)=z1E*dN1_dcsi(cs,et)+z2E*dN2_dcsi(cs,et)+z3E*dN3_dcsi(cs,et)+z4E*dN4_dcsi(cs,et);
    dz_deta(:,l)=z1E*dN1_deta(cs,et)+z2E*dN2_deta(cs,et)+z3E*dN3_deta(cs,et)+z4E*dN4_deta(cs,et);
end
detJE=sqrt((dy_dcsi.*dz_deta-dy_deta.*dz_dcsi).^2+(-dx_dcsi.*dz_deta+dx_deta.*dz_dcsi).^2+(dx_dcsi.*dy_deta-dx_deta.*dy_dcsi).^2);
M_e=zeros(12,12,size(detJE,1));
I=zeros(78*size(detJE,1),1);
J=zeros(78*size(detJE,1),1);
m_ij=J;
el_ij=m_ij;
BTB=zeros(12,12,size(detJE,1),4);
LEN_tot=0;
for l=1:size(detJE,1)
    for m=1:4
        B=[N1pg(l,m)*eye(3),N2pg(l,m)*eye(3),N3pg(l,m)*eye(3),N4pg(l,m)*eye(3)];
        BTB(:,:,l,m)=B.'*B;
    end
    M_e(:,:,l)=detJE(l,1)*BTB(:,:,l,1)+detJE(l,2)*BTB(:,:,l,2)+detJE(l,3)*BTB(:,:,l,3)+detJE(l,4)*BTB(:,:,l,4);
    [IND_I,IND_J,Mv]=find(tril(M_e(:,:,l)));
    LEN=length(IND_I);
    I(LEN_tot+(1:LEN)',1)=reshape(ENGINE_FACES_DOFS(l,IND_I),[],1);
    J(LEN_tot+(1:LEN)',1)=reshape(ENGINE_FACES_DOFS(l,IND_J),[],1);
    m_ij(LEN_tot+(1:LEN)',1)=Mv(:);
    el_ij(LEN_tot+(1:LEN)',1)=l;
    LEN_tot=LEN_tot+LEN;
end
I=I(1:LEN_tot);
J=J(1:LEN_tot);
m_ij=m_ij(1:LEN_tot);
el_ij=el_ij(1:LEN_tot);
i1=I;
j1=J;
i1(I<J)=J(I<J);
j1(I<J)=I(I<J);
I=i1;
J=j1;
ME=sparse(I,J,m_ij,size(DOFs_MAP,1),size(DOFs_MAP,1));ME=ME+triu(ME.',1);ME=(ME+ME.')/2;
[In,Jn]=find(diag(diag(ME)));
ME=ME+(speye(size(ME))+sparse(In,Jn,-1*ones(size(In)),size(ME,1),size(ME,1)));
%Projection operators RBF
%DZ->E
DISTANCE_EDZ=zeros(size(interface_E_coord,1),size(interface_DZ_coord,1));
for n=1:size(interface_DZ_coord,1)
    DISTANCE_EDZ(:,n)=sqrt(diag((interface_E_coord-ones(size(interface_E_coord,1),1)*interface_DZ_coord(n,:))*(interface_E_coord-ones(size(interface_E_coord,1),1)*interface_DZ_coord(n,:)).'));
end
DISTANCE_EE=zeros(size(interface_E_coord,1));
for n=1:size(interface_E_coord,1)
    DISTANCE_EE(:,n)=sqrt(diag((interface_E_coord-ones(size(interface_E_coord,1),1)*interface_E_coord(n,:))*(interface_E_coord-ones(size(interface_E_coord,1),1)*interface_E_coord(n,:)).'));
end
DISTANCE_DZDZ=zeros(size(interface_DZ_coord,1),size(interface_DZ_coord,1));
for n=1:size(interface_DZ_coord,1)
    DISTANCE_DZDZ(:,n)=sqrt(diag((interface_DZ_coord-ones(size(interface_DZ_coord,1),1)*interface_DZ_coord(n,:))*(interface_DZ_coord-ones(size(interface_DZ_coord,1),1)*interface_DZ_coord(n,:)).'));
end
%reference Radii evaluation
connectivity_matrix_E=zeros(9,size(interface_E_coord,1));
for nod=1:size(interface_E_coord,1)
    [Iel,~]=find(nod==FACE_intE(:,2:5));
    [neighbourhood_el]=(FACE_intE(Iel,2:5));
    neighbourhood_el=unique(neighbourhood_el(:));
    neighbourhood_el=setdiff(neighbourhood_el,nod);
    connectivity_matrix_E(1:(length(neighbourhood_el)+1),nod)=[nod;neighbourhood_el];
    for l=1:9
        if connectivity_matrix_E(l,nod)~=0
            Neig_dist_E(l,nod)=DISTANCE_EE(nod,connectivity_matrix_E(l,nod));
        else
            Neig_dist_E(l,nod)=0;
        end
    end
end

ELintidDZ=FACE_intdz;
for k=1:size(FACE_intdz,1)
    for l=1:size(FACE_intdz,2)-1
        ELintidDZ(k,l+1)=find(interface_DZ==FACE_intdz(k,l+1));
    end
end
connectivity_matrix_DZ=zeros(9,size(interface_DZ_coord,1));
for nod=1:size(interface_DZ_coord,1)
    [Iel,~]=find((nod)==ELintidDZ(:,2:5));
    [neighbourhood_el]=(ELintidDZ(Iel,2:5));
    neighbourhood_el=unique(neighbourhood_el(:));
    neighbourhood_el=setdiff(neighbourhood_el,(nod));
    connectivity_matrix_DZ(1:(length(neighbourhood_el)+1),nod)=[(nod);neighbourhood_el];
    for l=1:9
        if connectivity_matrix_DZ(l,nod)~=0
            Neig_dist_DZ(l,nod)=DISTANCE_DZDZ(nod,connectivity_matrix_DZ(l,nod));
        else
            Neig_dist_DZ(l,nod)=0;
        end
    end
end
radii_DZ=max(Neig_dist_DZ)/sqrt(2);
radii_E=max(Neig_dist_E)/sqrt(2);
radii_DZDZ=repmat(radii_DZ,size(DISTANCE_DZDZ,1),1);
PHI_MM_DZ=(1-DISTANCE_DZDZ./radii_DZDZ).^4.*(1-DISTANCE_DZDZ./radii_DZDZ>=0).*(1+4.*DISTANCE_DZDZ./radii_DZDZ);
radii_EDZ=repmat(radii_DZ,size(DISTANCE_EE,1),1);
PHI_NM_DZ=(1-DISTANCE_EDZ./radii_EDZ).^4.*(1-DISTANCE_EDZ./radii_EDZ>=0).*(1+4.*DISTANCE_EDZ./radii_EDZ);
radii_EE=repmat(radii_E,size(DISTANCE_EE,1),1);
PHI_MM_E=(1-DISTANCE_EE./radii_EE).^4.*(1-DISTANCE_EE./radii_EE>=0).*(1+4.*DISTANCE_EE./radii_EE);
radii_DZE=repmat(radii_E,size(DISTANCE_EDZ,2),1);
PHI_NM_E=(1-DISTANCE_EDZ.'./radii_DZE).^4.*(1-DISTANCE_EDZ.'./radii_DZE>=0).*(1+4.*DISTANCE_EDZ.'./radii_DZE);
Pi_g_DZ=PHI_NM_DZ*(PHI_MM_DZ\ones(size(PHI_MM_DZ,1),1));
Pi_g_DZ=repmat(Pi_g_DZ,1,size(PHI_MM_DZ,1));
Pi_g_DZ(Pi_g_DZ==0)=1;
Pi_g_E=PHI_NM_E/PHI_MM_E*ones(size(PHI_MM_E,1),1);
Pi_g_E=repmat(Pi_g_E,1,size(PHI_MM_E,1));
Pi_g_E(Pi_g_E==0)=1;
Pr_node_EDZ=(PHI_NM_DZ/PHI_MM_DZ)./(Pi_g_DZ);
Pr_node_DZE=(PHI_NM_E/PHI_MM_E)./(Pi_g_E);
PrDZE=zeros(size(MDZ,1),size(ME,1));
for nn=1:size(Pr_node_DZE,1);
    for mm=1:size(Pr_node_DZE,2)
        PrDZE(3*(interface_DZ(nn)-1)+(1:3),3*(mm-1)+(1:3)) = Pr_node_DZE(nn,mm)*eye(3);
    end
end
PrEDZ=zeros(size(MDZ,1),size(ME,1)).';
for nn=1:size(Pr_node_DZE,1);
    for mm=1:size(Pr_node_DZE,2)
        PrEDZ(3*(mm-1)+(1:3),3*(interface_DZ(nn)-1)+(1:3)) = Pr_node_EDZ(mm,nn)*eye(3);
    end
end
PrEDZ(end-sum(DOFs_MAP(:,1)==0)+1:end,end-sum(DOFs_MAP(:,1)==0)+1:end)=eye(sum(DOFs_MAP(:,1)==0));
PrDZE(end-sum(DOFs_MAP(:,1)==0)+1:end,end-sum(DOFs_MAP(:,1)==0)+1:end)=eye(sum(DOFs_MAP(:,1)==0));
Intergrid_switch=0;
if Intergrid_switch==1
    Q_E=MDZ*PrDZE/ME;
else
    Q_E=PrEDZ';
end
% K_def_compliance=PrEDZ'*K_interface*PrEDZ;

K_def=sparse(Q_E*K_interface*PrEDZ);
LOAD_MATRIX_def=sparse(Q_E*LOAD_MATRIX);
Recovery_Matrix_def=sparse(Recovery_Matrix*PrEDZ);
FAN_Center_Recovery_vector_def=sparse(FAN_Center_Recovery_vector*PrEDZ);
x1=COORD(DESIGN_ZONE_ELEMENT_NODES(:,1),1);
y1=COORD(DESIGN_ZONE_ELEMENT_NODES(:,1),2);
z1=COORD(DESIGN_ZONE_ELEMENT_NODES(:,1),3);
x2=COORD(DESIGN_ZONE_ELEMENT_NODES(:,2),1);
y2=COORD(DESIGN_ZONE_ELEMENT_NODES(:,2),2);
z2=COORD(DESIGN_ZONE_ELEMENT_NODES(:,2),3);
x3=COORD(DESIGN_ZONE_ELEMENT_NODES(:,3),1);
y3=COORD(DESIGN_ZONE_ELEMENT_NODES(:,3),2);
z3=COORD(DESIGN_ZONE_ELEMENT_NODES(:,3),3);
x4=COORD(DESIGN_ZONE_ELEMENT_NODES(:,4),1);
y4=COORD(DESIGN_ZONE_ELEMENT_NODES(:,4),2);
z4=COORD(DESIGN_ZONE_ELEMENT_NODES(:,4),3);
x5=COORD(DESIGN_ZONE_ELEMENT_NODES(:,5),1);
y5=COORD(DESIGN_ZONE_ELEMENT_NODES(:,5),2);
z5=COORD(DESIGN_ZONE_ELEMENT_NODES(:,5),3);
x6=COORD(DESIGN_ZONE_ELEMENT_NODES(:,6),1);
y6=COORD(DESIGN_ZONE_ELEMENT_NODES(:,6),2);
z6=COORD(DESIGN_ZONE_ELEMENT_NODES(:,6),3);
x7=COORD(DESIGN_ZONE_ELEMENT_NODES(:,7),1);
y7=COORD(DESIGN_ZONE_ELEMENT_NODES(:,7),2);
z7=COORD(DESIGN_ZONE_ELEMENT_NODES(:,7),3);
x8=COORD(DESIGN_ZONE_ELEMENT_NODES(:,8),1);
y8=COORD(DESIGN_ZONE_ELEMENT_NODES(:,8),2);
z8=COORD(DESIGN_ZONE_ELEMENT_NODES(:,8),3);
%stiffness matrix
N1=@(csi,eta,mu) 1/8*(1-csi).*(1-eta).*(1-mu);
N2=@(csi,eta,mu) 1/8*(1+csi).*(1-eta).*(1-mu);
N3=@(csi,eta,mu) 1/8*(1+csi).*(1+eta).*(1-mu);
N4=@(csi,eta,mu) 1/8*(1-csi).*(1+eta).*(1-mu);
N5=@(csi,eta,mu) 1/8*(1-csi).*(1-eta).*(1+mu);
N6=@(csi,eta,mu) 1/8*(1+csi).*(1-eta).*(1+mu);
N7=@(csi,eta,mu) 1/8*(1+csi).*(1+eta).*(1+mu);
N8=@(csi,eta,mu) 1/8*(1-csi).*(1+eta).*(1+mu);
E0=210000;
Emin=0.01;
nu=0.29;
%matrerial constant for unitary young modul
C1=(1-nu)/(1+nu)/(1-2*nu); %E1/E
C2=nu/(1+nu)/(1-2*nu); %E2/E
C3=1/2/(1+nu); %G/E
D=[C1 C2 C2 0 0 0
    C2 C1 C2 0 0 0
    C2 C2 C1 0 0 0
    0 0 0 C3 0 0
    0 0 0 0 C3 0
    0 0 0 0 0 C3];
dN1_dcsi=@(csi,eta,mu) -1/8*(1-eta).*(1-mu);
dN1_deta=@(csi,eta,mu) -1/8*(1-csi).*(1-mu);
dN1_dmu=@(csi,eta,mu) -1/8*(1-csi).*(1-eta);
%
dN2_dcsi=@(csi,eta,mu) 1/8*(1-eta).*(1-mu);
dN2_deta=@(csi,eta,mu) -1/8*(1+csi).*(1-mu);
dN2_dmu=@(csi,eta,mu) -1/8*(1+csi).*(1-eta);
%
dN3_dcsi=@(csi,eta,mu) 1/8*(1+eta).*(1-mu);
dN3_deta=@(csi,eta,mu) 1/8*(1+csi).*(1-mu);
dN3_dmu=@(csi,eta,mu) -1/8*(1+csi).*(1+eta);
%
dN4_dcsi=@(csi,eta,mu) -1/8*(1+eta).*(1-mu);
dN4_deta=@(csi,eta,mu) 1/8*(1-csi).*(1-mu);
dN4_dmu=@(csi,eta,mu) -1/8*(1-csi).*(1+eta);
%
dN5_dcsi=@(csi,eta,mu) -1/8*(1-eta).*(1+mu);
dN5_deta=@(csi,eta,mu) -1/8*(1-csi).*(1+mu);
dN5_dmu=@(csi,eta,mu) +1/8*(1-csi).*(1-eta);
%
dN6_dcsi=@(csi,eta,mu) 1/8*(1-eta).*(1+mu);
dN6_deta=@(csi,eta,mu) -1/8*(1+csi).*(1+mu);
dN6_dmu=@(csi,eta,mu) +1/8*(1+csi).*(1-eta);
%
dN7_dcsi=@(csi,eta,mu) 1/8*(1+eta).*(1+mu);
dN7_deta=@(csi,eta,mu) 1/8*(1+csi).*(1+mu);
dN7_dmu=@(csi,eta,mu) +1/8*(1+csi).*(1+eta);
%
dN8_dcsi=@(csi,eta,mu) -1/8*(1+eta).*(1+mu);
dN8_deta=@(csi,eta,mu) 1/8*(1-csi).*(1+mu);
dN8_dmu=@(csi,eta,mu) +1/8*(1-csi).*(1+eta);
Csi=[-1/sqrt(3),+1/sqrt(3),1/sqrt(3),-1/sqrt(3),-1/sqrt(3),+1/sqrt(3),1/sqrt(3),-1/sqrt(3)];
Eta=[-1/sqrt(3),-1/sqrt(3),1/sqrt(3),1/sqrt(3),-1/sqrt(3),-1/sqrt(3),1/sqrt(3),1/sqrt(3)];
Mu=[-1/sqrt(3),-1/sqrt(3),-1/sqrt(3),-1/sqrt(3),1/sqrt(3),1/sqrt(3),1/sqrt(3),1/sqrt(3)];
dx_dcsi=zeros(length(x1),8);
dx_deta=dx_dcsi;
dx_dmu=dx_dcsi;
dy_dcsi=dx_dcsi;
dy_deta=dx_dcsi;
dy_dmu=dx_dcsi;
dz_dcsi=dx_dcsi;
dz_deta=dx_dcsi;
dz_dmu=dx_dcsi;
for l=1:8
    
    cs=Csi(l);
    et=Eta(l);
    mui=Mu(l);
    dN1_dCsi(:,l)=dN1_dcsi(cs,et,mui);
    dN2_dCsi(:,l)=dN2_dcsi(cs,et,mui);
    dN3_dCsi(:,l)=dN3_dcsi(cs,et,mui);
    dN4_dCsi(:,l)=dN4_dcsi(cs,et,mui);
    dN5_dCsi(:,l)=dN5_dcsi(cs,et,mui);
    dN6_dCsi(:,l)=dN6_dcsi(cs,et,mui);
    dN7_dCsi(:,l)=dN7_dcsi(cs,et,mui);
    dN8_dCsi(:,l)=dN8_dcsi(cs,et,mui);
    dN1_dEta(:,l)=dN1_deta(cs,et,mui);
    dN2_dEta(:,l)=dN2_deta(cs,et,mui);
    dN3_dEta(:,l)=dN3_deta(cs,et,mui);
    dN4_dEta(:,l)=dN4_deta(cs,et,mui);
    dN5_dEta(:,l)=dN5_deta(cs,et,mui);
    dN6_dEta(:,l)=dN6_deta(cs,et,mui);
    dN7_dEta(:,l)=dN7_deta(cs,et,mui);
    dN8_dEta(:,l)=dN8_deta(cs,et,mui);
    dN1_dMu(:,l)=dN1_dmu(cs,et,mui);
    dN2_dMu(:,l)=dN2_dmu(cs,et,mui);
    dN3_dMu(:,l)=dN3_dmu(cs,et,mui);
    dN4_dMu(:,l)=dN4_dmu(cs,et,mui);
    dN5_dMu(:,l)=dN5_dmu(cs,et,mui);
    dN6_dMu(:,l)=dN6_dmu(cs,et,mui);
    dN7_dMu(:,l)=dN7_dmu(cs,et,mui);
    dN8_dMu(:,l)=dN8_dmu(cs,et,mui);
    dx_dcsi(:,l)=x1*dN1_dcsi(cs,et,mui)+x2*dN2_dcsi(cs,et,mui)+x3*dN3_dcsi(cs,et,mui)+x4*dN4_dcsi(cs,et,mui)+x5*dN5_dcsi(cs,et,mui)+x6*dN6_dcsi(cs,et,mui)+x7*dN7_dcsi(cs,et,mui)+x8*dN8_dcsi(cs,et,mui);
    dx_deta(:,l)=x1*dN1_deta(cs,et,mui)+x2*dN2_deta(cs,et,mui)+x3*dN3_deta(cs,et,mui)+x4*dN4_deta(cs,et,mui)+x5*dN5_deta(cs,et,mui)+x6*dN6_deta(cs,et,mui)+x7*dN7_deta(cs,et,mui)+x8*dN8_deta(cs,et,mui);
    dx_dmu(:,l)=x1*dN1_dmu(cs,et,mui)+x2*dN2_dmu(cs,et,mui)+x3*dN3_dmu(cs,et,mui)+x4*dN4_dmu(cs,et,mui)+x5*dN5_dmu(cs,et,mui)+x6*dN6_dmu(cs,et,mui)+x7*dN7_dmu(cs,et,mui)+x8*dN8_dmu(cs,et,mui);
    dy_dcsi(:,l)=y1*dN1_dcsi(cs,et,mui)+y2*dN2_dcsi(cs,et,mui)+y3*dN3_dcsi(cs,et,mui)+y4*dN4_dcsi(cs,et,mui)+y5*dN5_dcsi(cs,et,mui)+y6*dN6_dcsi(cs,et,mui)+y7*dN7_dcsi(cs,et,mui)+y8*dN8_dcsi(cs,et,mui);
    dy_deta(:,l)=y1*dN1_deta(cs,et,mui)+y2*dN2_deta(cs,et,mui)+y3*dN3_deta(cs,et,mui)+y4*dN4_deta(cs,et,mui)+y5*dN5_deta(cs,et,mui)+y6*dN6_deta(cs,et,mui)+y7*dN7_deta(cs,et,mui)+y8*dN8_deta(cs,et,mui);
    dy_dmu(:,l)=y1*dN1_dmu(cs,et,mui)+y2*dN2_dmu(cs,et,mui)+y3*dN3_dmu(cs,et,mui)+y4*dN4_dmu(cs,et,mui)+y5*dN5_dmu(cs,et,mui)+y6*dN6_dmu(cs,et,mui)+y7*dN7_dmu(cs,et,mui)+y8*dN8_dmu(cs,et,mui);
    dz_dcsi(:,l)=z1*dN1_dcsi(cs,et,mui)+z2*dN2_dcsi(cs,et,mui)+z3*dN3_dcsi(cs,et,mui)+z4*dN4_dcsi(cs,et,mui)+z5*dN5_dcsi(cs,et,mui)+z6*dN6_dcsi(cs,et,mui)+z7*dN7_dcsi(cs,et,mui)+z8*dN8_dcsi(cs,et,mui);
    dz_deta(:,l)=z1*dN1_deta(cs,et,mui)+z2*dN2_deta(cs,et,mui)+z3*dN3_deta(cs,et,mui)+z4*dN4_deta(cs,et,mui)+z5*dN5_deta(cs,et,mui)+z6*dN6_deta(cs,et,mui)+z7*dN7_deta(cs,et,mui)+z8*dN8_deta(cs,et,mui);
    dz_dmu(:,l)=z1*dN1_dmu(cs,et,mui)+z2*dN2_dmu(cs,et,mui)+z3*dN3_dmu(cs,et,mui)+z4*dN4_dmu(cs,et,mui)+z5*dN5_dmu(cs,et,mui)+z6*dN6_dmu(cs,et,mui)+z7*dN7_dmu(cs,et,mui)+z8*dN8_dmu(cs,et,mui);
end
A11=dy_deta.*dz_dmu-dy_dmu.*dz_deta;
A22=dz_dmu.*dx_dcsi-dx_dmu.*dz_dcsi;
A33=dx_dcsi.*dy_deta-dx_deta.*dy_dcsi;
A12=dx_dmu.*dz_deta-dx_deta.*dz_dmu;
A23=dx_dmu.*dy_dcsi-dx_dcsi.*dy_dmu;
A31=dy_dcsi.*dz_deta-dz_dcsi.*dy_deta;
A21=dz_dcsi.*dy_dmu-dy_dcsi.*dz_dmu;
A32=dz_dcsi.*dx_deta-dz_deta.*dx_dcsi;
A13=dx_deta.*dy_dmu-dy_deta.*dx_dmu;
detJ=dx_dcsi.*A11+dy_dcsi.*A12+dz_dcsi.*A13;
A2=sum(detJ,2)*2;
for l=1:8
    dN1_dx(:,l)=dN1_dCsi(l)*A11(:,l)./detJ(:,l)+dN1_dEta(l)*A21(:,l)./detJ(:,l)+dN1_dMu(l)*A31(:,l)./detJ(:,l);
    dN2_dx(:,l)=dN2_dCsi(l)*A11(:,l)./detJ(:,l)+dN2_dEta(l)*A21(:,l)./detJ(:,l)+dN2_dMu(l)*A31(:,l)./detJ(:,l);
    dN3_dx(:,l)=dN3_dCsi(l)*A11(:,l)./detJ(:,l)+dN3_dEta(l)*A21(:,l)./detJ(:,l)+dN3_dMu(l)*A31(:,l)./detJ(:,l);
    dN4_dx(:,l)=dN4_dCsi(l)*A11(:,l)./detJ(:,l)+dN4_dEta(l)*A21(:,l)./detJ(:,l)+dN4_dMu(l)*A31(:,l)./detJ(:,l);
    dN5_dx(:,l)=dN5_dCsi(l)*A11(:,l)./detJ(:,l)+dN5_dEta(l)*A21(:,l)./detJ(:,l)+dN5_dMu(l)*A31(:,l)./detJ(:,l);
    dN6_dx(:,l)=dN6_dCsi(l)*A11(:,l)./detJ(:,l)+dN6_dEta(l)*A21(:,l)./detJ(:,l)+dN6_dMu(l)*A31(:,l)./detJ(:,l);
    dN7_dx(:,l)=dN7_dCsi(l)*A11(:,l)./detJ(:,l)+dN7_dEta(l)*A21(:,l)./detJ(:,l)+dN7_dMu(l)*A31(:,l)./detJ(:,l);
    dN8_dx(:,l)=dN8_dCsi(l)*A11(:,l)./detJ(:,l)+dN8_dEta(l)*A21(:,l)./detJ(:,l)+dN8_dMu(l)*A31(:,l)./detJ(:,l);
    dN1_dy(:,l)=dN1_dCsi(l)*A12(:,l)./detJ(:,l)+dN1_dEta(l)*A22(:,l)./detJ(:,l)+dN1_dMu(l)*A32(:,l)./detJ(:,l);
    dN2_dy(:,l)=dN2_dCsi(l)*A12(:,l)./detJ(:,l)+dN2_dEta(l)*A22(:,l)./detJ(:,l)+dN2_dMu(l)*A32(:,l)./detJ(:,l);
    dN3_dy(:,l)=dN3_dCsi(l)*A12(:,l)./detJ(:,l)+dN3_dEta(l)*A22(:,l)./detJ(:,l)+dN3_dMu(l)*A32(:,l)./detJ(:,l);
    dN4_dy(:,l)=dN4_dCsi(l)*A12(:,l)./detJ(:,l)+dN4_dEta(l)*A22(:,l)./detJ(:,l)+dN4_dMu(l)*A32(:,l)./detJ(:,l);
    dN5_dy(:,l)=dN5_dCsi(l)*A12(:,l)./detJ(:,l)+dN5_dEta(l)*A22(:,l)./detJ(:,l)+dN5_dMu(l)*A32(:,l)./detJ(:,l);
    dN6_dy(:,l)=dN6_dCsi(l)*A12(:,l)./detJ(:,l)+dN6_dEta(l)*A22(:,l)./detJ(:,l)+dN6_dMu(l)*A32(:,l)./detJ(:,l);
    dN7_dy(:,l)=dN7_dCsi(l)*A12(:,l)./detJ(:,l)+dN7_dEta(l)*A22(:,l)./detJ(:,l)+dN7_dMu(l)*A32(:,l)./detJ(:,l);
    dN8_dy(:,l)=dN8_dCsi(l)*A12(:,l)./detJ(:,l)+dN8_dEta(l)*A22(:,l)./detJ(:,l)+dN8_dMu(l)*A32(:,l)./detJ(:,l);
    dN1_dz(:,l)=dN1_dCsi(l)*A13(:,l)./detJ(:,l)+dN1_dEta(l)*A23(:,l)./detJ(:,l)+dN1_dMu(l)*A33(:,l)./detJ(:,l);
    dN2_dz(:,l)=dN2_dCsi(l)*A13(:,l)./detJ(:,l)+dN2_dEta(l)*A23(:,l)./detJ(:,l)+dN2_dMu(l)*A33(:,l)./detJ(:,l);
    dN3_dz(:,l)=dN3_dCsi(l)*A13(:,l)./detJ(:,l)+dN3_dEta(l)*A23(:,l)./detJ(:,l)+dN3_dMu(l)*A33(:,l)./detJ(:,l);
    dN4_dz(:,l)=dN4_dCsi(l)*A13(:,l)./detJ(:,l)+dN4_dEta(l)*A23(:,l)./detJ(:,l)+dN4_dMu(l)*A33(:,l)./detJ(:,l);
    dN5_dz(:,l)=dN5_dCsi(l)*A13(:,l)./detJ(:,l)+dN5_dEta(l)*A23(:,l)./detJ(:,l)+dN5_dMu(l)*A33(:,l)./detJ(:,l);
    dN6_dz(:,l)=dN6_dCsi(l)*A13(:,l)./detJ(:,l)+dN6_dEta(l)*A23(:,l)./detJ(:,l)+dN6_dMu(l)*A33(:,l)./detJ(:,l);
    dN7_dz(:,l)=dN7_dCsi(l)*A13(:,l)./detJ(:,l)+dN7_dEta(l)*A23(:,l)./detJ(:,l)+dN7_dMu(l)*A33(:,l)./detJ(:,l);
    dN8_dz(:,l)=dN8_dCsi(l)*A13(:,l)./detJ(:,l)+dN8_dEta(l)*A23(:,l)./detJ(:,l)+dN8_dMu(l)*A33(:,l)./detJ(:,l);
end
Ke_E=zeros(24,24,size(detJ,1));
I=zeros(300*size(detJ,1),1);
J=zeros(300*size(detJ,1),1);
k_ij=J;
el_ij=k_ij;
BTDB=zeros(24,24,size(detJ,1),8);
for l=1:size(detJ,1)
    for m=1:8
        B=[diag([dN1_dx(l,m);dN1_dy(l,m);dN1_dz(l,m)]) ,diag([dN2_dx(l,m);dN2_dy(l,m);dN2_dz(l,m)]), diag([dN3_dx(l,m);dN3_dy(l,m);dN3_dz(l,m)]),diag([dN4_dx(l,m);dN4_dy(l,m);dN4_dz(l,m)]),...
            diag([dN5_dx(l,m);dN5_dy(l,m);dN5_dz(l,m)]), diag([dN6_dx(l,m);dN6_dy(l,m);dN6_dz(l,m)]), diag([dN7_dx(l,m);dN7_dy(l,m);dN7_dz(l,m)]), diag([dN8_dx(l,m);dN8_dy(l,m);dN8_dz(l,m)])
            [dN1_dy(l,m) dN1_dx(l,m) 0; 0 dN1_dz(l,m) dN1_dy(l,m);dN1_dz(l,m) 0 dN1_dx(l,m)], [dN2_dy(l,m) dN2_dx(l,m) 0; 0 dN2_dz(l,m) dN2_dy(l,m);dN2_dz(l,m) 0 dN2_dx(l,m)]...
            ,[dN3_dy(l,m) dN3_dx(l,m) 0; 0 dN3_dz(l,m) dN3_dy(l,m);dN3_dz(l,m) 0 dN3_dx(l,m)],[dN4_dy(l,m) dN4_dx(l,m) 0; 0 dN4_dz(l,m) dN4_dy(l,m);dN4_dz(l,m) 0 dN4_dx(l,m)],...
            [dN5_dy(l,m) dN5_dx(l,m) 0; 0 dN5_dz(l,m) dN5_dy(l,m);dN5_dz(l,m) 0 dN5_dx(l,m)],[dN6_dy(l,m) dN6_dx(l,m) 0; 0 dN6_dz(l,m) dN6_dy(l,m);dN6_dz(l,m) 0 dN6_dx(l,m)],...
            [dN7_dy(l,m) dN7_dx(l,m) 0; 0 dN7_dz(l,m) dN7_dy(l,m);dN7_dz(l,m) 0 dN7_dx(l,m)],[dN8_dy(l,m) dN8_dx(l,m) 0; 0 dN8_dz(l,m) dN8_dy(l,m);dN8_dz(l,m) 0 dN8_dx(l,m)],];
        BTDB(:,:,l,m)=B.'*D*B;
    end
    Ke_E(:,:,l)=detJ(l,1)*BTDB(:,:,l,1)+detJ(l,2)*BTDB(:,:,l,2)+detJ(l,3)*BTDB(:,:,l,3)+detJ(l,4)*BTDB(:,:,l,4)+detJ(l,5)*BTDB(:,:,l,5)+detJ(l,6)*BTDB(:,:,l,6)+detJ(l,7)*BTDB(:,:,l,7)+detJ(l,8)*BTDB(:,:,l,8);
    [IND_I,IND_J,Kv]=find(tril(Ke_E(:,:,l)));
    I(300*(l-1)+(1:300)',1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,IND_I),[],1);
    J(300*(l-1)+(1:300)',1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,IND_J),[],1);
    k_ij(300*(l-1)+(1:300)',1)=Kv(:);
    el_ij(300*(l-1)+(1:300)',1)=l;
end
i1=I;
j1=J;
i1(I<J)=J(I<J);
j1(I<J)=I(I<J);
I=i1;
J=j1;
% BCs
% fixed DOFs in the region x>=4000&Z==max(Z)
% BCs
% fixed DOFs in the region x>=4000&Z==max(Z)
fixed_nodes=DESIGN_ZONE_COORD(DESIGN_ZONE_COORD(:,2)>=4000&DESIGN_ZONE_COORD(:,4)>=max(DESIGN_ZONE_COORD(:,4))-10,1);
fixed_dofs=[3*(fixed_nodes-1)+1;3*(fixed_nodes-1)+2;3*(fixed_nodes-1)+3];
fixed_dofs=sort(fixed_dofs);
All_DOFs=(1:length(unique(DESIGN_ZONE_ELEMENT_DOFS(:)))).';
free_dofs=setdiff(All_DOFs,fixed_dofs);
Lagrangian_multip=new_DOFs_MAP(end-sum(new_DOFs_MAP(:,1)==0)+1:end,4);
disp_DOF=sort(setdiff(free_dofs,Lagrangian_multip));
free_dofs=union(free_dofs,Lagrangian_multip);
All_DOFs=union(All_DOFs,Lagrangian_multip);
U=zeros(length(All_DOFs),1);
% Ax_dz=(-x1dz-x2dz+x3dz+x4dz)/4;
% Bx_dz=(-x1dz+x2dz+x3dz-x4dz)/4;
% Cx_dz=(x1dz-x2dz+x3dz-x4dz)/4;
% Dx_dz=(x1dz+x2dz+x3dz+x4dz)/4;
% Ay_dz=(-y1dz-y2dz+y3dz+y4dz)/4;
% By_dz=(-y1dz+y2dz+y3dz-y4dz)/4;
% Cy_dz=(y1dz-y2dz+y3dz-y4dz)/4;
% Dy_dz=(y1dz+y2dz+y3dz+y4dz)/4;
% Az_dz=(-z1dz-z2dz+z3dz+z4dz)/4;
% Bz_dz=(-z1dz+z2dz+z3dz-z4dz)/4;
% Cz_dz=(z1dz-z2dz+z3dz-z4dz)/4;
% Dz_dz=(z1dz+z2dz+z3dz+z4dz)/4;
% csi_eta=zeros(length(Ax_dz),size(interface_E_coord,1));
% options = optimoptions('fminunc','Algorithm','trust-region','GradObj','off');
% for nn=1:length(Ax_dz)
%     for mm=1:size(interface_E_coord,1)
%
%         distance_squared=@(x) (Ax_dz(nn)*x(1)+Bx_dz(nn)*x(2)+Cx_dz(nn)*x(1).*x(2)+Dx_dz(nn)-interface_E_coord(mm,1)).^2+(Ay_dz(nn)*x(1)+By_dz(nn)*x(2)+Cy_dz(nn)*x(1).*x(2)+Dy_dz(nn)-interface_E_coord(mm,2)).^2+(Az_dz(nn)*x(1)+Bz_dz(nn)*x(2)+Cz_dz(nn)*x(1).*x(2)+Dz_dz(nn)-interface_E_coord(mm,3)).^2;
%         grad=@(x) [2*(Ax_dz(nn)*x(1)+Bx_dz(nn)*x(2)+Cx_dz(nn)*x(1).*x(2)+Dx_dz(nn)-interface_E_coord(mm,1))*(Ax_dz(nn)+Cx_dz(nn)*x(2))+2*(Ay_dz(nn)*x(1)+By_dz(nn)*x(2)+Cy_dz(nn)*x(1).*x(2)+Dy_dz(nn)-interface_E_coord(mm,2))*(Ay_dz(nn)+Cy_dz(nn)*x(2)+2*(Az_dz(nn)*x(1)+Bz_dz(nn)*x(2)+Cz_dz(nn)*x(1).*x(2)+Dz_dz(nn)-interface_E_coord(mm,3))*(Az_dz(nn)+Cz_dz(nn)*x(2)))
%             2*(Ax_dz(nn)*x(1)+Bx_dz(nn)*x(2)+Cx_dz(nn)*x(1).*x(2)+Dx_dz(nn)-interface_E_coord(mm,1))*(Bx_dz(nn)+Cx_dz(nn)*x(1))+2*(Ay_dz(nn)*x(1)+By_dz(nn)*x(2)+Cy_dz(nn)*x(1).*x(2)+Dy_dz(nn)-interface_E_coord(mm,2))*(By_dz(nn)+Cy_dz(nn)*x(1)+2*(Az_dz(nn)*x(1)+Bz_dz(nn)*x(2)+Cz_dz(nn)*x(1).*x(2)+Dz_dz(nn)-interface_E_coord(mm,3))*(Bz_dz(nn)+Cz_dz(nn)*x(1)))];
%
%         csi_eta((nn-1)*2+(1:2),mm)=fminunc(@(x) distance_squared(x),[0;0],options);
%         disp(['node =',num2str(mm),' Element =',num2str(nn)])
%     end
% end
%identify planar element that are easier
% planar_element=find(~(Cx_dz|Cy_dz|Cz_dz));
% Identification_mat=zeros(3*length(Ax_dz));
% Identification_vector=zeros(3*length(Ax_dz),1);
% for nn=1:length(Ax_dz)
%     Identification_mat(3*(nn-1)+(1:3),3*(nn-1)+(1:3))=[Ax_dz(nn) Bx_dz(nn) Cx_dz(nn)
%         Ay_dz(nn) By_dz(nn) Cy_dz(nn)
%         Az_dz(nn) Bz_dz(nn) Cz_dz(nn)];
%     Identification_vector(3*(nn-1)+(1:3),1)=[Dx_dz(nn);Dy_dz(nn);Dz_dz(nn)];
% end
% for mm=1:size(interface_E_coord,1)
%     pool_matrix(:,mm)=repmat(interface_E_coord(mm,:).',length(Ax_dz),1)-Identification_vector;
% end
% Identification_mat=Identification_mat+eye(size(Identification_mat))*0.000000000001;
% csi_eta_csieta=pinv(Identification_mat)*pool_matrix;
%identify master candidates

% ELintidDZ=FACE_intdz;
% for k=1:size(FACE_intdz,1)
%     for l=1:size(FACE_intdz,2)-1
%         ELintidDZ(k,l+1)=find(interface_DZ==FACE_intdz(k,l+1));
%     end
% end
% hold on; h24_patch = patch('Vertices',interface_DZ_coord,'Faces',ELintidDZ(:,2:5),'FaceVertexCData',[1 1 1],'FaceColor','flat'); axis equal; axis off; drawnow;hold off
%
% hold on
% scatter3(interface_E_coord(:,1),interface_E_coord(:,2),interface_E_coord(:,3),'fill')
% for mm=1:size(intercoord,1)
%     text(interface_E_coord(mm,1),interface_E_coord(mm,2),interface_E_coord(mm,3),num2str(mm));
% end
% for mm=1:size(interface_DZ_coord,1)
%     text(interface_DZ_coord(mm,1),interface_DZ_coord(mm,2),interface_DZ_coord(mm,3),num2str(mm));
% end
%nearest_point

% ELintid=ELint;
% for k=1:size(ELint,1)
%     for l=1:size(ELint,2)-1
%        ELintid(k,l+1)=find(intercoord(:,1)==ELint(k,l+1));
%     end
% end
% hold on; h24_patch = patch('Vertices',intercoord(:,2:end),'Faces',ELintid(:,2:5),'FaceVertexCData',[1 1 1],'FaceColor','flat'); axis equal; axis off; drawnow;hold off
% for mm=1:size(intercoord,1)
%     text(intercoord(mm,2),intercoord(mm,3),intercoord(mm,4),num2str(intercoord(mm,1)));
% end
% hold on
% scatter3(interface_DZ_coord(:,1),interface_DZ_coord(:,2),interface_DZ_coord(:,3),'fill')
% hold on; h24_patch = patch('Vertices',COORD(:,[1,2,3]),'Faces',ELEMENT,'FaceVertexCData',0.3*[1 1 1],'FaceColor','flat'); axis equal; axis off; drawnow;hold off
% for mm=1:size(COORD,1)
%     text(COORD(mm,1),COORD(mm,2),COORD(mm,3),num2str(mm));
% end
%projection_matrix
%% INITIALIZE ITERATION
volfrac=0.3;
penal=3; ft=1;
x = volfrac*ones(size(DESIGN_ZONE_ELEMENT_NODES,1),1);
xPhys = x(:);

% xPhys(x1>=3000&x4<=5000&z1<=min(z1)+dz)=0;
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
maxoutit  = 10000;
kkttol  = 1e-3;

%
%%%% The iterations start:
kktnorm = kkttol+1000;
outit = 0;
%% Preparing Filter
rmin=2;
center_node_coordinate=[(x1+x2+x3+x4+x5+x6+x7+x8)/8,(y1+y2+y3+y4+y5+y6+y7+y8)/8,(z1+z2+z3+z4+z5+z6+z7+z8)/8];
kit=0;
iH=zeros(length(xPhys)^2,1);
jH=zeros(length(xPhys)^2,1);
sH=zeros(length(xPhys)^2,1);
for k=1:length(xPhys)
    for j=1:length(xPhys)
        kit=kit+1;
        iH(kit)=k;
        jH(kit)=j;
        sH(kit) = max(0,rmin*80-sqrt((center_node_coordinate(k,1)-center_node_coordinate(j,1))^2+(center_node_coordinate(k,2)-center_node_coordinate(j,2))^2+(center_node_coordinate(k,3)-center_node_coordinate(j,3))^2));
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% FE-ANALYSIS
E=Emin+xPhys(:).^penal*(E0-Emin);
E=E(el_ij);
dK=penal*(E0-Emin)*xPhys(:).^(penal-1);
dK=dK(el_ij);
Lind=(1:length(xPhys(:))).';
Lind=Lind(el_ij);
dK=k_ij.*dK;
sK=k_ij.*E;
tic
K = sparse(I,J,sK,size(K_def,1),size(K_def,2)); K=K+triu(K.',1); K=(K+K.')/2;
K=K_def+K;
F=[LOAD_MATRIX_def(:,1)];
tic
Res=K(free_dofs,free_dofs)\F(free_dofs,:);
inversion_time=toc
U(free_dofs)=Res(:,1);
% Res=(K(free_dofs,free_dofs)).'\F(free_dofs,:);
load('dKDXs');
tic
dK_dxU=zeros(length(U),length(xPhys(:)));
xPhys=xPhys(:);
for xdim=1:length(xPhys(:))
    dK_dxU(:,xdim)=(penal*(E0-Emin)*xPhys(xdim).^(penal-1))*dK_DXs(xdim).mat*U;
end
dk_dxu_eval=toc
%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
Utip=Recovery_Matrix_def*U;
perfo=0;
DTSFC_DU=0;
for s=1:20
    tip_vector=Gamma(s).mat*(U_static(:,1)+Utip);
    R=rms(tip_vector);
    perfo=perfo+Gamma(s).lambda*R;
    N=length(tip_vector);
    R=ones(N,1)*R;
    DTSFC_DU=DTSFC_DU+Gamma(s).lambda/N*Recovery_Matrix_def.'*Gamma(s).mat.'*(tip_vector./R);
end
tic
RES2=K(free_dofs,free_dofs).'\[DTSFC_DU(free_dofs),FAN_Center_Recovery_vector_def(free_dofs).'];
psi_tild=zeros(size(K,1),2);
adjoint_evaluation=toc
psi_tild(free_dofs,:)=RES2;
FAN_Ax_disp=U_FAN_static(1)+FAN_Center_Recovery_vector_def*U;
dperfo=-psi_tild(:,1).'*dK_dxU;
dfandisp=-psi_tild(:,2).'*dK_dxU;
dv = A2/2;
% FILTERING/MODIFICATION OF SENSITIVITIES
if ft == 1
    dperfo(:) = H*(x(:).*dperfo.')./Hs./max(1e-3,x(:));
    dfandisp(:)=H*(x(:).*dfandisp.')./Hs./max(1e-3,x(:));
elseif ft == 2
    dperfo = H*(dperfo.'./Hs);
    dfandisp=H*(dfandisp.'./Hs);
    dv(:) = H*(dv(:)./Hs);
end
f0val=10*perfo;
fval=[(sum(xPhys(:).*A2/2) - volfrac*sum(A2/2))/sum(xPhys(:).*A2/2);1*(-0.38-FAN_Ax_disp)/0.38];%
df0dx=10*reshape(dperfo,[],1);
dfdx=[dv(:)'/sum(xPhys(:).*A2/2);-1*dfandisp(:)'/0.38];%
innerit=0;
outvector1 = [outeriter innerit f0val fval'];
outvector2 = xval;
Fac_color=zeros(6*length(xval),1);
for nn=1:length(xval)
    Fac_color(6*(nn-1)+(1:6))=(1-xval(nn))*ones(6,1);
end
%color filter
FACES_shown=FACES(Fac_color<=0.9,:);
Fac_color=Fac_color(Fac_color<=0.9);
h = figure(1); set(h,'Color',[1 1 1]);
clf
hold on; h24_patch = patch('Vertices',COORD,'Faces',FACES_shown(:,2:end),'FaceVertexCData',Fac_color*[1 1 1],'FaceColor','flat'); axis equal; axis off;
quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
text(1000,0,0,'x')
text(0,1000,0,'y')
text(0,0,1000,'z')
view([27.6 18]);
drawnow;
hold off
fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',outeriter,perfo, ...
    sum(xPhys(:).*A2/2)/sum(A2/2),kktnorm);
%
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
    %     xPhys(x1>=2500&x4<=4500&z1==min(z1))=0;
    %% FE-ANALYSIS
    E=Emin+xPhys(:).^penal*(E0-Emin);
    E=E(el_ij);
    dK=penal*(E0-Emin)*xPhys(:).^(penal-1);
    dK=dK(el_ij);
    Lind=(1:length(xPhys(:))).';
    Lind=Lind(el_ij);
    dK=k_ij.*dK;
    sK=k_ij.*E;
    K = sparse(I,J,sK,size(K_def,1),size(K_def,2)); K=K+triu(K.',1); K=(K+K.')/2;
    K=K_def+K;
    F=[LOAD_MATRIX_def(:,1)];
    Res=K(free_dofs,free_dofs)\F(free_dofs,:);
    U(free_dofs)=Res(:,1);
    %     Res=(K(free_dofs,free_dofs)).'\F(free_dofs,:);
    %     [Res,flag,relres,iter,resvec]= gmres(K(free_dofs,free_dofs),F(free_dofs,:),[],1e-4,1000);
    %     flag
    %     iter
    %     relres
    U(free_dofs)=Res(:,1);
    dK_dxU=zeros(length(U),length(xPhys(:)));
    xPhys=xPhys(:);
    for xdim=1:length(xPhys(:))
        dK_dxU(:,xdim)=(penal*(E0-Emin)*xPhys(xdim).^(penal-1))*dK_DXs(xdim).mat*U;
    end
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    Utip=Recovery_Matrix_def*U;
    perfo=0;
    DTSFC_DU=0;
    for s=1:20
        tip_vector=Gamma(s).mat*(U_static(:,1)+Utip);
        R=rms(tip_vector);
        perfo=perfo+Gamma(s).lambda*R;
        N=length(tip_vector);
        R=ones(N,1)*R;
        DTSFC_DU=DTSFC_DU+Gamma(s).lambda/N*Recovery_Matrix_def.'*Gamma(s).mat.'*(tip_vector./R);
    end
    RES2=K(free_dofs,free_dofs).'\[DTSFC_DU(free_dofs),FAN_Center_Recovery_vector_def(free_dofs).'];
    psi_tild=zeros(size(K,1),2);
    psi_tild(free_dofs,:)=RES2;
    FAN_Ax_disp=U_FAN_static(1)+FAN_Center_Recovery_vector_def*U;
    dperfo=-psi_tild(:,1).'*dK_dxU;
    dfandisp=-psi_tild(:,2).'*dK_dxU;
    dv = A2/2;
    % FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dperfo(:) = H*(x(:).*dperfo.')./Hs./max(1e-3,x(:));
        dfandisp(:)=H*(x(:).*dfandisp.')./Hs./max(1e-3,x(:));
    elseif ft == 2
        dperfo = H*(dperfo.'./Hs);
        dfandisp=H*(dfandisp.'./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
    f0val=10*perfo;
    fval=[(sum(xPhys(:).*A2/2) - volfrac*sum(A2/2))/sum(xPhys(:).*A2/2);1*(-0.38-FAN_Ax_disp)/0.38];%
    df0dx=10*reshape(dperfo,[],1);
    dfdx=[dv(:)'/sum(xPhys(:).*A2/2);-1*dfandisp(:)'/0.38];%
    %% %% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
    outvector1 = [outeriter innerit f0val fval(:)'];
    outvector2 = xval';
    xPhys=xval';
    Fac_color=zeros(6*length(xval),1);
    for nn=1:length(xval)
        Fac_color(6*(nn-1)+(1:6))=(1-xval(nn))*ones(6,1);
    end
    %color filter
    FACES_shown=FACES(Fac_color<=0.3,:);
    Fac_color=Fac_color(Fac_color<=0.3);
    h = figure(1); set(h,'Color',[1 1 1]);
    clf
    hold on; h24_patch = patch('Vertices',COORD,'Faces',FACES_shown(:,2:end),'FaceVertexCData',Fac_color*[1 1 1],'FaceColor','flat'); axis equal; axis off;
    quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
    quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
    quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
    text(1000,0,0,'x')
    text(0,1000,0,'y')
    text(0,0,1000,'z')
    view([27.6 18]);
    drawnow;
    hold off
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',outeriter,perfo, ...
        sum(xPhys(:).*A2/2)/sum(A2/2),kktnorm);
    %
    figure(2)
    hold on
    scatter(outeriter,f0val,'fill','k')
    figure(3)
    hold on
    scatter(outeriter,FAN_Ax_disp,'fill','b')
end
% Average dimension
for nn=1:size(disp_DOF,1)
    Grid_motion(new_DOFs_MAP(disp_DOF(nn),1),new_DOFs_MAP(disp_DOF(nn),3))=U(disp_DOF(nn));
end
MeanSize = mean(max(COORD)-min(COORD));
% Maximum motion
[MaxGridMotion,I] = max(abs(real(Grid_motion(:))));
% MaxGridMotion = MaxGridMotion*sign(real(Grid_motion(I)));
% New grid location
%     NewGridPositionResidual = COORD(:,2:4) + 0.5*MeanSize*real(GridMotion_residual)/MaxGridMotionResidual;
% Color displacement
ColorOdsR = sqrt(sum(real(Grid_motion)'.^2))';
ColorOdsR = ColorOdsR - min(ColorOdsR); ColorOdsR = ColorOdsR/max(ColorOdsR);
% Plot

Final_coord=zeros(size(DESIGN_ZONE_COORD));
for k=1:size(DESIGN_ZONE_COORD,1)
    Final_coord(k,:)=DESIGN_ZONE_COORD(k,:)+0.5*MeanSize*[0,U(3*(k-1)+1),U(3*(k-1)+2),U(3*k)]/MaxGridMotion;
end
figure
hold on; h25_patch = patch('Vertices',Final_coord(:,2:end),'Faces',FACES_shown(:,2:end),'CData',ColorOdsR,'FaceColor','interp','SpecularColorReflectance',0.3);
colormap jet;
quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
text(1000,0,0,'x')
text(0,1000,0,'y')
text(0,0,1000,'z')
view([27.6 18]);

