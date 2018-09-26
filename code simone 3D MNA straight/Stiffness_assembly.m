function [I,J,k_ij,Is,Js,DBs]=Stiffness_assembly(COORD,DESIGN_ZONE_ELEMENT_NODES,DESIGN_ZONE_ELEMENT_DOFS)
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
%
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
Is=zeros(120*size(detJ,1),1);
Js=zeros(120*size(detJ,1),1);
DBs=zeros(120*size(detJ,1),1);
k_ij=J;
el_ij=k_ij;
BTDB=zeros(24,24,size(detJ,1),8);
DBe=zeros(6,24,size(detJ,1),8);
DB=zeros(6,24,size(detJ,1));
for l=1:size(detJ,1)
    for m=1:8
        B=[diag([dN1_dx(l,m);dN1_dy(l,m);dN1_dz(l,m)]) ,diag([dN2_dx(l,m);dN2_dy(l,m);dN2_dz(l,m)]), diag([dN3_dx(l,m);dN3_dy(l,m);dN3_dz(l,m)]),diag([dN4_dx(l,m);dN4_dy(l,m);dN4_dz(l,m)]),...
            diag([dN5_dx(l,m);dN5_dy(l,m);dN5_dz(l,m)]), diag([dN6_dx(l,m);dN6_dy(l,m);dN6_dz(l,m)]), diag([dN7_dx(l,m);dN7_dy(l,m);dN7_dz(l,m)]), diag([dN8_dx(l,m);dN8_dy(l,m);dN8_dz(l,m)])
            [dN1_dy(l,m) dN1_dx(l,m) 0; 0 dN1_dz(l,m) dN1_dy(l,m);dN1_dz(l,m) 0 dN1_dx(l,m)], [dN2_dy(l,m) dN2_dx(l,m) 0; 0 dN2_dz(l,m) dN2_dy(l,m);dN2_dz(l,m) 0 dN2_dx(l,m)]...
            ,[dN3_dy(l,m) dN3_dx(l,m) 0; 0 dN3_dz(l,m) dN3_dy(l,m);dN3_dz(l,m) 0 dN3_dx(l,m)],[dN4_dy(l,m) dN4_dx(l,m) 0; 0 dN4_dz(l,m) dN4_dy(l,m);dN4_dz(l,m) 0 dN4_dx(l,m)],...
            [dN5_dy(l,m) dN5_dx(l,m) 0; 0 dN5_dz(l,m) dN5_dy(l,m);dN5_dz(l,m) 0 dN5_dx(l,m)],[dN6_dy(l,m) dN6_dx(l,m) 0; 0 dN6_dz(l,m) dN6_dy(l,m);dN6_dz(l,m) 0 dN6_dx(l,m)],...
            [dN7_dy(l,m) dN7_dx(l,m) 0; 0 dN7_dz(l,m) dN7_dy(l,m);dN7_dz(l,m) 0 dN7_dx(l,m)],[dN8_dy(l,m) dN8_dx(l,m) 0; 0 dN8_dz(l,m) dN8_dy(l,m);dN8_dz(l,m) 0 dN8_dx(l,m)],];
        DBe(:,:,l,m)=D*B;
        BTDB(:,:,l,m)=B.'*DBe(:,:,l,m);
    end
    Ke_E(:,:,l)=detJ(l,1)*BTDB(:,:,l,1)+detJ(l,2)*BTDB(:,:,l,2)+detJ(l,3)*BTDB(:,:,l,3)+detJ(l,4)*BTDB(:,:,l,4)+detJ(l,5)*BTDB(:,:,l,5)+detJ(l,6)*BTDB(:,:,l,6)+detJ(l,7)*BTDB(:,:,l,7)+detJ(l,8)*BTDB(:,:,l,8);
    [IND_I,IND_J,Kv]=find(tril(Ke_E(:,:,l)));
    I(300*(l-1)+(1:300)',1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,IND_I),[],1);
    J(300*(l-1)+(1:300)',1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,IND_J),[],1);
    k_ij(300*(l-1)+(1:300)',1)=E0*Kv(:);
    el_ij(300*(l-1)+(1:300)',1)=l;
    DB(:,:,l)=detJ(l,1)*DBe(:,:,l,1)+detJ(l,2)*DBe(:,:,l,2)+detJ(l,3)*DBe(:,:,l,3)+detJ(l,4)*DBe(:,:,l,4)+detJ(l,5)*DBe(:,:,l,5)+detJ(l,6)*DBe(:,:,l,6)+detJ(l,7)*DBe(:,:,l,7)+detJ(l,8)*DBe(:,:,l,8);
    DB(:,:,l)= DB(:,:,l)/sum(detJ(l,:));
    [IND_si,IND_sj,dbs]=find(DB(:,:,l));
    Is(120*(l-1)+(1:120),1)=reshape(6*(l-1)+IND_si,[],1);
    Js(120*(l-1)+(1:120),1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,IND_sj),[],1);
    DBs(120*(l-1)+(1:120))=E0*dbs(:);
end
i1=I;
j1=J;
i1(I<J)=J(I<J);
j1(I<J)=I(I<J);
I=i1;
J=j1;