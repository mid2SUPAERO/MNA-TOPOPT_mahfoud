function [FEM_structure]= Stiffness_components3D_compliance(FEM_structure,E0)
%% Extract inputs
ELEMENT=FEM_structure.ELEMENT;
COORD=FEM_structure.COORD;
ELEMENT=[(1:size(ELEMENT,1))',ELEMENT];
ELEMENTo=FEM_structure.Old_element;
ELEMENTo=[(1:size(ELEMENTo,1))',ELEMENTo];
COORDo=FEM_structure.Old_coord;
new_DOFs_MAP=FEM_structure.new_DOFs_MAP;

%% Evaluate unit stiffness matrix
first_node_index=ELEMENTo(:,2);
second_node_index=ELEMENTo(:,3);
third_node_index=ELEMENTo(:,4);
fourth_node_index=ELEMENTo(:,5);
fivth_node_index=ELEMENTo(:,6);
sixth_node_index=ELEMENTo(:,7);
seventh_node_index=ELEMENTo(:,8);
eigth_node_index=ELEMENTo(:,9);
DESIGN_ZONE_ELEMENT_NODES=[first_node_index,second_node_index,third_node_index,fourth_node_index,fivth_node_index,sixth_node_index,seventh_node_index,eigth_node_index];
x1o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,1),1);
y1o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,1),2);
z1o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,1),3);
x2o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,2),1);
y2o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,2),2);
z2o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,2),3);
x3o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,3),1);
y3o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,3),2);
z3o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,3),3);
x4o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,4),1);
y4o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,4),2);
z4o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,4),3);
x5o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,5),1);
y5o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,5),2);
z5o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,5),3);
x6o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,6),1);
y6o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,6),2);
z6o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,6),3);
x7o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,7),1);
y7o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,7),2);
z7o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,7),3);
x8o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,8),1);
y8o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,8),2);
z8o=COORDo(DESIGN_ZONE_ELEMENT_NODES(:,8),3);
DESIGN_ZONE_COORD=[(1:size(COORD,1)).',COORD];
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
% E0=210000;
% Emin=E0/1000;
rho=8e-9;
Emin=E0/1e6;
Emin0=Emin;
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
Ngpn=zeros(48,48);
% L=zeros(3,24,8);
for l=1:8
    cs=Csi(l);
    et=Eta(l);
    mui=Mu(l);
    Ngpn(6*(l-1)+(1:6),:)=[N1(cs,et,mui)*eye(6),N2(cs,et,mui)*eye(6),N3(cs,et,mui)*eye(6),N4(cs,et,mui)*eye(6),N5(cs,et,mui)*eye(6),N6(cs,et,mui)*eye(6),N7(cs,et,mui)*eye(6),N8(cs,et,mui)*eye(6)];
%     L(:,:,l)=sparse([N1(cs,et,mui)*eye(3),N2(cs,et,mui)*eye(3),N3(cs,et,mui)*eye(3),N4(cs,et,mui)*eye(3),N5(cs,et,mui)*eye(3),N6(cs,et,mui)*eye(3),N7(cs,et,mui)*eye(3),N8(cs,et,mui)*eye(3)]);
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
%     LTL(:,:,l)=L(:,:,l)'*L(:,:,l);
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
[Indgp,indn,Nn]=find(Ngpn);
numb_comp_nn=length(Nn);
In=zeros(size(detJ,1)*numb_comp_nn,1);
node_id_n=[repmat(DESIGN_ZONE_ELEMENT_NODES(:,1),1,6),repmat(DESIGN_ZONE_ELEMENT_NODES(:,2),1,6),repmat(DESIGN_ZONE_ELEMENT_NODES(:,3),1,6),repmat(DESIGN_ZONE_ELEMENT_NODES(:,4),1,6),repmat(DESIGN_ZONE_ELEMENT_NODES(:,5),1,6),repmat(DESIGN_ZONE_ELEMENT_NODES(:,6),1,6),repmat(DESIGN_ZONE_ELEMENT_NODES(:,7),1,6),repmat(DESIGN_ZONE_ELEMENT_NODES(:,8),1,6)];
node_id_s=[repmat(1:6,size(DESIGN_ZONE_ELEMENT_NODES,1),8)];
tl=0;
Jn=In;
Nnn=In;
I=zeros(300*size(detJ,1),1);
J=zeros(300*size(detJ,1),1);
Is=zeros(120*size(detJ,1),1);
Js=zeros(120*size(detJ,1),1);
DBs=zeros(120*size(detJ,1),1);
% Im=zeros(108*size(detJ,1),1);
% Jm=zeros(108*size(detJ,1),1);
% M_ij=zeros(108*size(detJ,1),1);
% elm_ij=M_ij;
k_ij=J;
el_ij=k_ij;
BTDB=zeros(24,24,8);
IK0=zeros(24*24*size(detJ,1),1);
JK0=zeros(24*24*size(detJ,1),1);
KV0=IK0;
DB=zeros(6*8,24);
for l=1:size(detJ,1)
    for m=1:8
        
        B=[diag([dN1_dx(l,m);dN1_dy(l,m);dN1_dz(l,m)]) ,diag([dN2_dx(l,m);dN2_dy(l,m);dN2_dz(l,m)]), diag([dN3_dx(l,m);dN3_dy(l,m);dN3_dz(l,m)]),diag([dN4_dx(l,m);dN4_dy(l,m);dN4_dz(l,m)]),...
            diag([dN5_dx(l,m);dN5_dy(l,m);dN5_dz(l,m)]), diag([dN6_dx(l,m);dN6_dy(l,m);dN6_dz(l,m)]), diag([dN7_dx(l,m);dN7_dy(l,m);dN7_dz(l,m)]), diag([dN8_dx(l,m);dN8_dy(l,m);dN8_dz(l,m)])
            [dN1_dy(l,m) dN1_dx(l,m) 0; 0 dN1_dz(l,m) dN1_dy(l,m);dN1_dz(l,m) 0 dN1_dx(l,m)], [dN2_dy(l,m) dN2_dx(l,m) 0; 0 dN2_dz(l,m) dN2_dy(l,m);dN2_dz(l,m) 0 dN2_dx(l,m)]...
            ,[dN3_dy(l,m) dN3_dx(l,m) 0; 0 dN3_dz(l,m) dN3_dy(l,m);dN3_dz(l,m) 0 dN3_dx(l,m)],[dN4_dy(l,m) dN4_dx(l,m) 0; 0 dN4_dz(l,m) dN4_dy(l,m);dN4_dz(l,m) 0 dN4_dx(l,m)],...
            [dN5_dy(l,m) dN5_dx(l,m) 0; 0 dN5_dz(l,m) dN5_dy(l,m);dN5_dz(l,m) 0 dN5_dx(l,m)],[dN6_dy(l,m) dN6_dx(l,m) 0; 0 dN6_dz(l,m) dN6_dy(l,m);dN6_dz(l,m) 0 dN6_dx(l,m)],...
            [dN7_dy(l,m) dN7_dx(l,m) 0; 0 dN7_dz(l,m) dN7_dy(l,m);dN7_dz(l,m) 0 dN7_dx(l,m)],[dN8_dy(l,m) dN8_dx(l,m) 0; 0 dN8_dz(l,m) dN8_dy(l,m);dN8_dz(l,m) 0 dN8_dx(l,m)],];
        DBe(:,:,m)=D*B;
        BTDB(:,:,m)=B.'*DBe(:,:,m);
        DB(6*(m-1)+(1:6),:)=DBe(:,:,m);
    end
%     Me_E=rho*(detJ(l,1)*LTL(:,:,1)+detJ(l,2)*LTL(:,:,2)+detJ(l,3)*LTL(:,:,3)+detJ(l,4)*LTL(:,:,4)+detJ(l,5)*LTL(:,:,5)+detJ(l,6)*LTL(:,:,6)+detJ(l,7)*LTL(:,:,7)+detJ(l,8)*LTL(:,:,8));
    Ke_E=detJ(l,1)*BTDB(:,:,1)+detJ(l,2)*BTDB(:,:,2)+detJ(l,3)*BTDB(:,:,3)+detJ(l,4)*BTDB(:,:,4)+detJ(l,5)*BTDB(:,:,5)+detJ(l,6)*BTDB(:,:,6)+detJ(l,7)*BTDB(:,:,7)+detJ(l,8)*BTDB(:,:,8);
    [IND_I,IND_J,Kv]=find(tril(Ke_E));
%     [IND_IM,IND_JM,Mv]=find(tril(Me_E));
    I(300*(l-1)+(1:300)',1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,IND_I),[],1);
    J(300*(l-1)+(1:300)',1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,IND_J),[],1);
    k_ij(300*(l-1)+(1:300)',1)=Kv(:);
    el_ij(300*(l-1)+(1:300)',1)=l;
%     Im(108*(l-1)+(1:108)',1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,IND_IM),[],1);
%     Jm(108*(l-1)+(1:108)',1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,IND_JM),[],1);
%     M_ij(108*(l-1)+(1:108)',1)=Mv(:);
%     elm_ij(108*(l-1)+(1:108)',1)=l;
    [Ii,Ji,KV]=find((Ke_E));
%     IK0(64*(l-1)+(1:64)',1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,Ii),[],1);
    IK0(576*(l-1)+(1:576)',1)=Ii(:)+(24)*(l-1);
    JK0(576*(l-1)+(1:576)',1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,Ji),[],1);
    KV0(576*(l-1)+(1:576)',1)=KV(:);
    In(numb_comp_nn*(l-1)+(1:numb_comp_nn))=48*(l-1)+Indgp;
    Jn(numb_comp_nn*(l-1)+(1:numb_comp_nn))=6*(node_id_n(l,indn)-1)+node_id_s(l,indn);
    Nnn(numb_comp_nn*(l-1)+(1:numb_comp_nn))=Nn;
    [IND_si,IND_sj,dbs]=find(DB);
    Ncomp=length(IND_si);
    Is(tl+(1:Ncomp),1)=reshape(48*(l-1)+IND_si,[],1);
    Js(tl+(1:Ncomp),1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,IND_sj),[],1);
    DBs(tl+(1:Ncomp))=E0*dbs(:);
    tl=tl+Ncomp;
    disp(num2str(l/size(detJ,1)*100))
end
K0c=sparse((IK0),(JK0),(KV0),size(detJ,1)*24,max(max(DESIGN_ZONE_ELEMENT_DOFS)));
clear IK0 JK0 KV0

clear BTDB DBe
i1=I;
j1=J;
i1(I<J)=J(I<J);
j1(I<J)=I(I<J);
I=i1;
J=j1;
DBs=DBs(1:tl);
Is=Is(1:tl);
Js=Js(1:tl);
DB=sparse(Is,Js,DBs,48*size(DESIGN_ZONE_ELEMENT_DOFS,1),max(max(DESIGN_ZONE_ELEMENT_DOFS)));
Ngpn=sparse(In,Jn,Nnn,48*size(DESIGN_ZONE_ELEMENT_DOFS,1),6*max(max(DESIGN_ZONE_ELEMENT_NODES)));

% A faire le stress
% DB=sparse(Is,Js,DBs,12*size(DESIGN_ZONE_ELEMENT_DOFS,1),max(max(DESIGN_ZONE_ELEMENT_DOFS)));
% Ngpn=sparse(In,Jn,Nnn,12*size(DESIGN_ZONE_ELEMENT_DOFS,1),3*max(max(DESIGN_ZONE_ELEMENT_NODES)));
%% Apply BCs
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
ndof=max(max(DESIGN_ZONE_ELEMENT_DOFS));fixeddofs=fixed_dofs;
ND = ones(ndof,1); ND(fixeddofs) = 0; Null = (spdiags(ND,0,ndof,ndof));
psi_tild=zeros(ndof,3); Sol=zeros(length(FEM_structure.observation_dofs)+size(FEM_structure.K_interface,1),1); psi_tild_sol=[Sol,Sol,Sol];
free_dofs_engine_master=size(FEM_structure.K_interface,1)+find(ismember(FEM_structure.observation_dofs,free_dofs));
free_dofs_engine_master=[(1:size(FEM_structure.K_interface,1))';free_dofs_engine_master];
%% Gather  outputs
FEM_structure.k_ij=(k_ij);
FEM_structure.I=(I);
FEM_structure.J=(J);
% FEM_structure.DB=(DB);
% FEM_structure.Ngpn=(Ngpn);
FEM_structure.x1=x1;
FEM_structure.y1=y1;
FEM_structure.z1=z1;
FEM_structure.x2=x2;
FEM_structure.y2=y2;
FEM_structure.z2=z2;
FEM_structure.x3=x3;
FEM_structure.y3=y3;
FEM_structure.z3=z3;
FEM_structure.x4=x4;
FEM_structure.y4=y4;
FEM_structure.z4=z4;
FEM_structure.x5=x5;
FEM_structure.y5=y5;
FEM_structure.z5=z5;
FEM_structure.x6=x6;
FEM_structure.y6=y6;
FEM_structure.z6=z6;
FEM_structure.x7=x7;
FEM_structure.y7=y7;
FEM_structure.z7=z7;
FEM_structure.x8=x8;
FEM_structure.y8=y8;
FEM_structure.z8=z8;
FEM_structure.A2=A2;
FEM_structure.E0=E0;
FEM_structure.Emin=Emin;
FEM_structure.Emin0=Emin0;
FEM_structure.free_dofs=free_dofs;
FEM_structure.DESIGN_ZONE_ELEMENT_NODES=DESIGN_ZONE_ELEMENT_NODES;
FEM_structure.DESIGN_ZONE_ELEMENT_DOFS=DESIGN_ZONE_ELEMENT_DOFS;
FEM_structure.DESIGN_ZONE_COORD=DESIGN_ZONE_COORD;
FEM_structure.el_ij=el_ij;
FEM_structure.U=U;
FEM_structure.K0c=(K0c);
FEM_structure.x1o=x1o;
FEM_structure.y1o=y1o;
FEM_structure.z1o=z1o;
FEM_structure.x2o=x2o;
FEM_structure.y2o=y2o;
FEM_structure.z2o=z2o;
FEM_structure.x3o=x3o;
FEM_structure.y3o=y3o;
FEM_structure.z3o=z3o;
FEM_structure.x4o=x4o;
FEM_structure.y4o=y4o;
FEM_structure.z4o=z4o;
FEM_structure.x5o=x5o;
FEM_structure.y5o=y5o;
FEM_structure.z5o=z5o;
FEM_structure.x6o=x6o;
FEM_structure.y6o=y6o;
FEM_structure.z6o=z6o;
FEM_structure.x7o=x7o;
FEM_structure.y7o=y7o;
FEM_structure.z7o=z7o;
FEM_structure.x8o=x8o;
FEM_structure.y8o=y8o;
FEM_structure.z8o=z8o;
FEM_structure.Null=Null;
FEM_structure.psi_tild=psi_tild;
FEM_structure.free_dofs_engine_master=free_dofs_engine_master;
FEM_structure.Sol=Sol;
FEM_structure.psi_tild_sol=psi_tild_sol;
FEM_structure.DB=(DB);
FEM_structure.Ngpn=(Ngpn);
% FEM_structure.M_ij=(M_ij);
% FEM_structure.Im=(Im);
% FEM_structure.Jm=(Jm);
% FEM_structure.elm_ij=elm_ij;
