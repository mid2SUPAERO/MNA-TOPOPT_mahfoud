function [FEM_structure]= engine_DOFs_eliminationpar(FEM_structure)
%% get inputs
FACES=FEM_structure.FACES;
if size(FEM_structure.FACES,2)==4
    FACES=[(1:size(FACES,1))',FACES];
end
ELint=FEM_structure.ELint;
%% interpolation.
new_DOFs_MAP=zeros(3*size(FEM_structure.COORD,1)+sum(FEM_structure.DOFs_MAP(:,1)==0),4);
for k=1:size(FEM_structure.COORD,1)
    new_DOFs_MAP(3*k-2,:)=[k,k,1,3*k-2];
    new_DOFs_MAP(3*k-1,:)=[k,k,2,3*k-1];
    new_DOFs_MAP(3*k,:)=[k,k,3,3*k];
end
new_DOFs_MAP(end-sum(FEM_structure.DOFs_MAP(:,1)==0)+1:end,:)=[FEM_structure.DOFs_MAP(end-sum(FEM_structure.DOFs_MAP(:,1)==0)+1:end,1:3),(3*k+(1:sum(FEM_structure.DOFs_MAP(:,1)==0)))'];

interface_DZ=FEM_structure.NODE_SET.NID;
interface_DZ_coord=FEM_structure.COORD(interface_DZ,:);
[ID4]=ismember(FACES(:,[2:5]),interface_DZ);
IDF=ID4(:,1)&ID4(:,2)&ID4(:,3)&ID4(:,4);
FACE_intdz=FACES(IDF,:);
x1dz=FEM_structure.COORD(FACE_intdz(:,2),1);
y1dz=FEM_structure.COORD(FACE_intdz(:,2),2);
z1dz=FEM_structure.COORD(FACE_intdz(:,2),3);
x2dz=FEM_structure.COORD(FACE_intdz(:,3),1);
y2dz=FEM_structure.COORD(FACE_intdz(:,3),2);
z2dz=FEM_structure.COORD(FACE_intdz(:,3),3);
x3dz=FEM_structure.COORD(FACE_intdz(:,4),1);
y3dz=FEM_structure.COORD(FACE_intdz(:,4),2);
z3dz=FEM_structure.COORD(FACE_intdz(:,4),3);
x4dz=FEM_structure.COORD(FACE_intdz(:,5),1);
y4dz=FEM_structure.COORD(FACE_intdz(:,5),2);
z4dz=FEM_structure.COORD(FACE_intdz(:,5),3);
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
N1pg=zeros(length(x1dz),4);
N2pg=zeros(length(x1dz),4);
N3pg=zeros(length(x1dz),4);
N4pg=zeros(length(x1dz),4);
dN1_dCsi=zeros(length(x1dz),4);
dN2_dCsi=zeros(length(x1dz),4);
dN3_dCsi=zeros(length(x1dz),4);
dN4_dCsi=zeros(length(x1dz),4);
dN1_dEta=zeros(length(x1dz),4);
dN2_dEta=zeros(length(x1dz),4);
dN3_dEta=zeros(length(x1dz),4);
dN4_dEta=zeros(length(x1dz),4);
dx_dcsi=zeros(length(x1dz),4);
dx_deta=zeros(length(x1dz),4);
dy_dcsi=zeros(length(x1dz),4);
dy_deta=zeros(length(x1dz),4);
dz_dcsi=zeros(length(x1dz),4);
dz_deta=zeros(length(x1dz),4);

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
BTB=zeros(12,12,4);
LEN_tot=0;
for l=1:size(detJdz,1)
    for m=1:4
        B=[N1pg(l,m)*eye(3),N2pg(l,m)*eye(3),N3pg(l,m)*eye(3),N4pg(l,m)*eye(3)];
        BTB(:,:,m)=B.'*B;
    end
    M_e(:,:,l)=detJdz(l,1)*BTB(:,:,1)+detJdz(l,2)*BTB(:,:,2)+detJdz(l,3)*BTB(:,:,3)+detJdz(l,4)*BTB(:,:,4);
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
MDZ=MDZ+(speye(size(MDZ))+sparse(In,Jn,-1*ones(size(In)),3*size(FEM_structure.COORD,1)+sum(FEM_structure.DOFs_MAP(:,1)==0),3*size(FEM_structure.COORD,1)+sum(FEM_structure.DOFs_MAP(:,1)==0)));
%same for engine
interface_E=FEM_structure.intercoord(:,1);
interface_E_coord=FEM_structure.intercoord(:,2:4);
FACE_intE=ELint;
for k=1:size(ELint,1)
    for l=2:size(ELint,2)
        id_kl=ELint(k,l);
        FACE_intE(k,l)=find(interface_E==id_kl);
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
BTB=zeros(12,12,4);
LEN_tot=0;
for l=1:size(detJE,1)
    for m=1:4
        B=[N1pg(l,m)*eye(3),N2pg(l,m)*eye(3),N3pg(l,m)*eye(3),N4pg(l,m)*eye(3)];
        BTB(:,:,m)=B.'*B;
    end
    M_e(:,:,l)=detJE(l,1)*BTB(:,:,1)+detJE(l,2)*BTB(:,:,2)+detJE(l,3)*BTB(:,:,3)+detJE(l,4)*BTB(:,:,4);
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
ME=sparse(I,J,m_ij,size(FEM_structure.DOFs_MAP,1),size(FEM_structure.DOFs_MAP,1));ME=ME+triu(ME.',1);ME=(ME+ME.')/2;
[In,Jn]=find(diag(diag(ME)));
ME=ME+(speye(size(ME))+sparse(In,Jn,-1*ones(size(In)),size(ME,1),size(ME,1)));
%Projection operators RBF
%DZ->E
iE=1:size(interface_E_coord,1);
iD=1:size(interface_DZ_coord,1);
[IE,ID]=meshgrid(iE,iD);
IE=IE(:);ID=ID(:);
DISTANCE_EDZ=sqrt(sum((interface_E_coord(IE,:)-interface_DZ_coord(ID,:)).^2,2));
DISTANCE_EDZ=sparse(IE,ID,DISTANCE_EDZ,size(interface_E_coord,1),size(interface_DZ_coord,1));
[IE,ID]=meshgrid(iE,iE);
IE=IE(:);ID=ID(:);
DISTANCE_EE=sqrt(sum((interface_E_coord(IE,:)-interface_E_coord(ID,:)).^2,2));
DISTANCE_EE=sparse(IE,ID,DISTANCE_EE,size(interface_E_coord,1),size(interface_E_coord,1));
[IE,ID]=meshgrid(iD,iD);
IE=IE(:);ID=ID(:);
DISTANCE_DZDZ=sqrt(sum((interface_DZ_coord(IE,:)-interface_DZ_coord(ID,:)).^2,2));
DISTANCE_DZDZ=sparse(IE,ID,DISTANCE_DZDZ,size(interface_DZ_coord,1),size(interface_DZ_coord,1));
% DISTANCE_EE=zeros(size(interface_E_coord,1));
% for n=1:size(interface_E_coord,1)
%     DISTANCE_EE(:,n)=sqrt(diag((interface_E_coord-ones(size(interface_E_coord,1),1)*interface_E_coord(n,:))*(interface_E_coord-ones(size(interface_E_coord,1),1)*interface_E_coord(n,:)).'));
% end
% DISTANCE_DZDZ=zeros(size(interface_DZ_coord,1),size(interface_DZ_coord,1));
% for n=1:size(interface_DZ_coord,1)
%     DISTANCE_DZDZ(:,n)=sqrt(diag((interface_DZ_coord-ones(size(interface_DZ_coord,1),1)*interface_DZ_coord(n,:))*(interface_DZ_coord-ones(size(interface_DZ_coord,1),1)*interface_DZ_coord(n,:)).'));
% end
%reference Radii evaluation
connectivity_matrix_E=zeros(9,size(interface_E_coord,1));
Neig_dist_E=connectivity_matrix_E;
parfor nod=1:size(interface_E_coord,1)
    [Iel,~]=find(nod==FACE_intE(:,2:5));
    [neighbourhood_el]=(FACE_intE(Iel,2:5));
    neighbourhood_el=unique(neighbourhood_el(:));
    neighbourhood_el=setdiff(neighbourhood_el,nod);
    D_EE=DISTANCE_EE(nod,:);
    for nei=1:9
        if nei==1
%             connectivity_matrix_E(nei,nod)=[nod];
            Neig_dist_E(nei,nod)=D_EE(nod);
        elseif nei<=(length(neighbourhood_el)+1)
%             connectivity_matrix_E(nei,nod)=[neighbourhood_el(nei-1)];
            Neig_dist_E(nei,nod)=D_EE(neighbourhood_el(nei-1));
        end
    end
end

ELintidDZ=FACE_intdz;
for k=1:size(FACE_intdz,1)
    for l=2:size(FACE_intdz,2)
        ELintidDZ(k,l)=find(interface_DZ==FACE_intdz(k,l));
    end
end
connectivity_matrix_DZ=zeros(9,size(interface_DZ_coord,1));
Neig_dist_DZ=connectivity_matrix_DZ;
parfor nod=1:size(interface_DZ_coord,1)
    [Iel,~]=find((nod)==ELintidDZ(:,2:5));
    [neighbourhood_el]=(ELintidDZ(Iel,2:5));
    neighbourhood_el=unique(neighbourhood_el(:));
    neighbourhood_el=setdiff(neighbourhood_el,(nod));
    D_DZDZ=DISTANCE_DZDZ(nod,:);
    for nei=1:9
        if nei==1
%             connectivity_matrix_E(nei,nod)=[nod];
            Neig_dist_DZ(nei,nod)=D_DZDZ(nod);
        elseif nei<=(length(neighbourhood_el)+1)
%             connectivity_matrix_E(nei,nod)=[neighbourhood_el(nei-1)];
            Neig_dist_DZ(nei,nod)=D_DZDZ(neighbourhood_el(nei-1));
        end
    end
end
radii_DZ=max(Neig_dist_DZ)/sqrt(2);
radii_E=max(Neig_dist_E)/sqrt(2);
radii_DZDZ=repmat(radii_DZ,size(DISTANCE_DZDZ,1),1);
PHI_MM_DZ=sparse((1-DISTANCE_DZDZ./radii_DZDZ).^4.*(1-DISTANCE_DZDZ./radii_DZDZ>=0).*(1+4.*DISTANCE_DZDZ./radii_DZDZ));
radii_EDZ=repmat(radii_DZ,size(DISTANCE_EE,1),1);
PHI_NM_DZ=sparse((1-DISTANCE_EDZ./radii_EDZ).^4.*(1-DISTANCE_EDZ./radii_EDZ>=0).*(1+4.*DISTANCE_EDZ./radii_EDZ));
radii_EE=repmat(radii_E,size(DISTANCE_EE,1),1);
PHI_MM_E=sparse((1-DISTANCE_EE./radii_EE).^4.*(1-DISTANCE_EE./radii_EE>=0).*(1+4.*DISTANCE_EE./radii_EE));
radii_DZE=repmat(radii_E,size(DISTANCE_EDZ,2),1);
PHI_NM_E=sparse((1-DISTANCE_EDZ.'./radii_DZE).^4.*(1-DISTANCE_EDZ.'./radii_DZE>=0).*(1+4.*DISTANCE_EDZ.'./radii_DZE));
Pi_g_DZ=PHI_NM_DZ*(PHI_MM_DZ\ones(size(PHI_MM_DZ,1),1));
Pi_g_DZ=repmat(Pi_g_DZ,1,size(PHI_MM_DZ,1));
Pi_g_DZ(Pi_g_DZ==0)=1;
Pi_g_E=PHI_NM_E/PHI_MM_E*ones(size(PHI_MM_E,1),1);
Pi_g_E=repmat(Pi_g_E,1,size(PHI_MM_E,1));
Pi_g_E(Pi_g_E==0)=1;
Pr_node_EDZ=(PHI_NM_DZ/PHI_MM_DZ)./(Pi_g_DZ);
Pr_node_DZE=(PHI_NM_E/PHI_MM_E)./(Pi_g_E);
PrDZE=zeros(size(MDZ,1),size(ME,1));
for nn=1:size(Pr_node_DZE,1)
    for mm=1:size(Pr_node_DZE,2)
        PrDZE(3*(interface_DZ(nn)-1)+(1:3),3*(mm-1)+(1:3)) = Pr_node_DZE(nn,mm)*eye(3);
    end
end
PrEDZ=zeros(size(MDZ,1),size(ME,1)).';
for nn=1:size(Pr_node_DZE,1)
    for mm=1:size(Pr_node_DZE,2)
        PrEDZ(3*(mm-1)+(1:3),3*(interface_DZ(nn)-1)+(1:3)) = Pr_node_EDZ(mm,nn)*eye(3);
    end
end
PrEDZ(end-sum(FEM_structure.DOFs_MAP(:,1)==0)+1:end,end-sum(FEM_structure.DOFs_MAP(:,1)==0)+1:end)=eye(sum(FEM_structure.DOFs_MAP(:,1)==0));
PrDZE(end-sum(FEM_structure.DOFs_MAP(:,1)==0)+1:end,end-sum(FEM_structure.DOFs_MAP(:,1)==0)+1:end)=eye(sum(FEM_structure.DOFs_MAP(:,1)==0));
Intergrid_switch=0;
if Intergrid_switch==1
    Q_E=MDZ*PrDZE/ME;
else
    Q_E=PrEDZ';
end
% for each  DOF that is not expressed in function of engine DOF
% the elimination is not applie
%the final list of DOFs is [u_Ec;u_Do]
observation_dofs=find(~any(PrDZE,2));
coupling_dofs=find(any(PrDZE,2));
Pr_EM=sparse([coupling_dofs;observation_dofs],(1:(size(PrDZE,1)))',ones(size(PrDZE,1),1))*[sparse(PrDZE(coupling_dofs,:)),spalloc(length(coupling_dofs),length(observation_dofs),1);spalloc(length(observation_dofs),size(PrDZE,2),1),speye(length(observation_dofs))];
K_Engine=spalloc(size(PrDZE,2)+length(observation_dofs),size(PrDZE,2)+length(observation_dofs),size(FEM_structure.K_interface,1)^2);
K_Engine(1:size(PrDZE,2),1:size(PrDZE,2))=FEM_structure.K_interface;
Pr=sparse(PrEDZ);
Q_E=Pr';
K_def=sparse(Q_E*FEM_structure.K_interface*Pr);
LOAD_MATRIX_def=Q_E*FEM_structure.LOAD_MATRIX;LOAD_MATRIX_EM=[sparse(FEM_structure.LOAD_MATRIX);spalloc(length(observation_dofs),3,1)];
Recovery_Matrix_def=sparse(FEM_structure.Recovery_Matrix*Pr);Recovery_Matrix_EM=[sparse(FEM_structure.Recovery_Matrix),spalloc(size(FEM_structure.Recovery_Matrix,1),length(observation_dofs),1)];
FAN_Center_Recovery_vector_def=sparse(FEM_structure.FAN_Center_Recovery_vector*Pr);
RG=cell(20,1);
RG_EM=cell(20,1);
parfor s=1:20
    N=size(FEM_structure.Gamma(s).mat,1);
    RG{s}=sparse(FEM_structure.Gamma(s).lambda/N*Recovery_Matrix_def.'*FEM_structure.Gamma(s).mat.');
    RG_EM{s}=sparse(FEM_structure.Gamma(s).lambda/N*Recovery_Matrix_EM.'*FEM_structure.Gamma(s).mat.');
end
%% gather outputs
FEM_structure.K_def=K_def;
FEM_structure.LOAD_MATRIX_def=LOAD_MATRIX_def;
FEM_structure.Recovery_Matrix_def=Recovery_Matrix_def;
FEM_structure.FAN_Center_Recovery_vector_def=FAN_Center_Recovery_vector_def;
FEM_structure.new_DOFs_MAP=new_DOFs_MAP;
FEM_structure.RG=RG;
FEM_structure.Pr_EM=Pr_EM;
FEM_structure.K_Engine=K_Engine;
FEM_structure.observation_dofs=observation_dofs;
FEM_structure.LOAD_MATRIX_EM=LOAD_MATRIX_EM;
FEM_structure.Recovery_Matrix_EM=Recovery_Matrix_EM;
FEM_structure.RG_EM=RG_EM;

