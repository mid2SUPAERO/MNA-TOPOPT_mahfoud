clear
close all
load('substructure_results')
load('lin_pert10.dat.mat')
load('dzstr.mat')
load('connectivity_interface')
%the point at the same x have the same 1 and 3 diplacement + linear
%% interpolation.
%find the interpolation nodes of the Design ZONE
interface_DZ=NODE_SET.NID;
interface_DZ_coord=COORD(interface_DZ,:);
%nearest_point
X_disc_intercoord=unique(intercoord(:,2));
Master_nodes=zeros(size(X_disc_intercoord));
closest_DZnodes=zeros(size(X_disc_intercoord,1),2);
for l=1:length(X_disc_intercoord)
    Master_nodes_Z=intercoord(intercoord(:,2)==X_disc_intercoord(l),[1,4]);
    Master_nodes(l)=Master_nodes_Z(Master_nodes_Z(:,2)==max(Master_nodes_Z(:,2)),1);
    Slave_nodes(l,:)=reshape(setdiff(Master_nodes_Z(:,1),Master_nodes(l)),1,[]);
    distance=ones(size(interface_DZ_coord,1),1)*intercoord(intercoord(:,1)==Master_nodes(l),2:end)-interface_DZ_coord;
    distance=sqrt(distance(:,1).^2+distance(:,2).^2+distance(:,3).^2);
    [~,Id]=sort(distance);
    closest_DZnodes(l,:)=interface_DZ([Id(1),Id(2)]);
    weight(l,:)=[distance(Id(2))/(distance(Id(2))+distance(Id(1))),distance(Id(1))/(distance(Id(2))+distance(Id(1)))];
end
All_Slave=[Master_nodes,Slave_nodes];
Pr=zeros(size(K_interface,1),2*size(COORD,1)+sum(DOFs_MAP(:,1)==0));
for k=1:size(K_interface,1)
    Node_k=DOFs_MAP(k,1);
    if ~isempty(Node_k)
    [I_k,~]=find(All_Slave==Node_k);
    DOF_DZ=[2*(closest_DZnodes(I_k,:)-1)+1,2*(closest_DZnodes(I_k,:)-1)+2];
    Ww=[weight(I_k,:),weight(I_k,:)];
    Pr(k,DOF_DZ)=Ww;
    end
end
Pr(end-sum(DOFs_MAP(:,1)==0)+1:end,end-sum(DOFs_MAP(:,1)==0)+1:end)=eye(sum(DOFs_MAP(:,1)==0));
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
% hold on; h24_patch = patch('Vertices',COORD(:,[1,2,3]),'Faces',ELEMENT,'FaceVertexCData',0.3*[1 1 1],'FaceColor','flat'); axis equal; axis off; drawnow;hold off
% for mm=1:size(COORD,1)
%     text(COORD(mm,1),COORD(mm,2),COORD(mm,3),num2str(mm));
% end
%projection_matrix
new_DOFs_MAP=zeros(2*size(COORD,1)+sum(DOFs_MAP(:,1)==0),4);
for k=1:size(COORD,1)
new_DOFs_MAP(2*k-1,:)=[k,k,1,2*k-1];
new_DOFs_MAP(2*k,:)=[k,k,3,2*k];
end
new_DOFs_MAP(end-sum(DOFs_MAP(:,1)==0)+1:end,:)=[DOFs_MAP(end-sum(DOFs_MAP(:,1)==0)+1:end,1:3),(2*k+(1:sum(DOFs_MAP(:,1)==0)))'];
K_def=Pr.'*K_interface*Pr;
LOAD_MATRIX_def=Pr.'*LOAD_MATRIX;
Recovery_Matrix_def=Recovery_Matrix*Pr;
FAN_Center_Recovery_vector_def=FAN_Center_Recovery_vector*Pr;
%

DESIGN_ZONE_COORD=[(1:size(COORD,1)).',COORD];
El_id=(1:size(ELEMENT,1)).';
first_node_index=ELEMENT(:,1);
second_node_index=ELEMENT(:,2);
third_node_index=ELEMENT(:,3);
fourth_node_index=ELEMENT(:,4);
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
%stiffness matrix
N1=@(csi,eta) 1/4*(1-csi).*(1-eta);
N2=@(csi,eta) 1/4*(1-csi).*(1+eta);
N3=@(csi,eta) 1/4*(1+csi).*(1+eta);
N4=@(csi,eta) 1/4*(1+csi).*(1-eta);
E0=210000;
Emin=1;
nu=0.3;
%matrerial constant for unitary young modul
C1=(1-nu)/(1+nu)/(1-2*nu); %E1/E
C2=C1*nu/(1-nu); %E2/E
C3=1/2/(1+nu); %G/E
D=[C1 C2 0
    C2 C1 0
    0 0 C3];
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
    dN1_dCsi(:,l)=dN1_dcsi(cs,et);
    dN2_dCsi(:,l)=dN2_dcsi(cs,et);
    dN3_dCsi(:,l)=dN3_dcsi(cs,et);
    dN4_dCsi(:,l)=dN4_dcsi(cs,et);
    dN1_dEta(:,l)=dN1_deta(cs,et);
    dN2_dEta(:,l)=dN2_deta(cs,et);
    dN3_dEta(:,l)=dN3_deta(cs,et);
    dN4_dEta(:,l)=dN4_deta(cs,et);
    dx_dcsi(:,l)=x1*dN1_dcsi(cs,et)+x2*dN2_dcsi(cs,et)+x3*dN3_dcsi(cs,et)+x4*dN4_dcsi(cs,et);
    dx_deta(:,l)=x1*dN1_deta(cs,et)+x2*dN2_deta(cs,et)+x3*dN3_deta(cs,et)+x4*dN4_deta(cs,et);
    dz_dcsi(:,l)=z1*dN1_dcsi(cs,et)+z2*dN2_dcsi(cs,et)+z3*dN3_dcsi(cs,et)+z4*dN4_dcsi(cs,et);
    dz_deta(:,l)=z1*dN1_deta(cs,et)+z2*dN2_deta(cs,et)+z3*dN3_deta(cs,et)+z4*dN4_deta(cs,et);
    
end
detJ=dx_dcsi.*dz_deta-dx_deta.*dz_dcsi;
I11=dz_deta./detJ;
I12=-dz_dcsi./detJ;
I21=-dx_deta./detJ;
I22=dx_dcsi./detJ;
for l=1:4
    dN1_dx(:,l)=dN1_dCsi(l)*I11(:,l)+dN1_dEta(l)*I12(:,l);
    dN2_dx(:,l)=dN2_dCsi(l)*I11(:,l)+dN2_dEta(l)*I12(:,l);
    dN3_dx(:,l)=dN3_dCsi(l)*I11(:,l)+dN3_dEta(l)*I12(:,l);
    dN4_dx(:,l)=dN4_dCsi(l)*I11(:,l)+dN4_dEta(l)*I12(:,l);
    dN1_dz(:,l)=dN1_dCsi(l)*I21(:,l)+dN1_dEta(l)*I22(:,l);
    dN2_dz(:,l)=dN2_dCsi(l)*I21(:,l)+dN2_dEta(l)*I22(:,l);
    dN3_dz(:,l)=dN3_dCsi(l)*I21(:,l)+dN3_dEta(l)*I22(:,l);
    dN4_dz(:,l)=dN4_dCsi(l)*I21(:,l)+dN4_dEta(l)*I22(:,l);
    
end
Ke_E=zeros(8,8,size(detJ,1));
I=zeros(36*size(detJ,1),1);
J=zeros(36*size(detJ,1),1);
k_ij=J;
el_ij=k_ij;
for l=1:size(detJ,1)
    for m=1:4
        B=[dN1_dx(l,m) 0 dN2_dx(l,m) 0 dN3_dx(l,m) 0 dN4_dx(l,m) 0
            0 dN1_dz(l,m) 0 dN2_dz(l,m) 0 dN3_dz(l,m) 0 dN4_dz(l,m)
            dN1_dz(l,m) dN1_dx(l,m) dN2_dz(l,m) dN2_dx(l,m) dN3_dz(l,m) dN3_dx(l,m) dN4_dz(l,m) dN4_dx(l,m)];
        BTDB(:,:,l,m)=B.'*D*B;
    end
    Ke_E(:,:,l)=detJ(l,1)*BTDB(:,:,l,1)+detJ(l,2)*BTDB(:,:,l,2)+detJ(l,3)*BTDB(:,:,l,3)+detJ(l,4)*BTDB(:,:,l,4);
    [IND_I,IND_J,Kv]=find(tril(Ke_E(:,:,l)));
    I(36*(l-1)+(1:36)',1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,IND_I),[],1);
    J(36*(l-1)+(1:36)',1)=reshape(DESIGN_ZONE_ELEMENT_DOFS(l,IND_J),[],1);
    k_ij(36*(l-1)+(1:36)',1)=Kv(:);
    el_ij(36*(l-1)+(1:36)',1)=l;
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
fixed_nodes=DESIGN_ZONE_COORD(DESIGN_ZONE_COORD(:,2)>=4000&DESIGN_ZONE_COORD(:,4)>=max(DESIGN_ZONE_COORD(:,4))-0.01,1);
fixed_dofs=[2*(fixed_nodes-1)+1;2*(fixed_nodes-1)+2];
fixed_dofs=sort(fixed_dofs);
All_DOFs=(1:length(unique(DESIGN_ZONE_ELEMENT_DOFS(:)))).';
free_dofs=setdiff(All_DOFs,fixed_dofs);
Lagrangian_multip=new_DOFs_MAP(end-sum(new_DOFs_MAP(:,1)==0)+1:end,4);
disp_DOF=free_dofs;
free_dofs=union(free_dofs,Lagrangian_multip);
All_DOFs=union(All_DOFs,Lagrangian_multip);
U=zeros(length(All_DOFs),1);


% %% Preparing Filter
% rmin=2;
% nelx=nx-1;
% nely=nz-1;
% iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
% jH = ones(size(iH));
% sH = zeros(size(iH));
% k = 0;
% center_node_coordinate=[(x1+x2+x3+x4)/4,(z1+z2+z3+z4)/4];
% for i1 = 1:nely
%     for j1 = 1:nelx
%         e1 = (i1-1)*nelx+j1;
%         for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nely)
%             for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nelx)
%                 e2 = (i2-1)*nelx+j2;
%                 k = k+1;
%                 iH(k) = e1;
%                 jH(k) = e2;
%                 sH(k) = max(0,rmin*mean(diff(X_dis))-sqrt((center_node_coordinate(e1,1)-center_node_coordinate(e2,1))^2+(center_node_coordinate(e1,2)-center_node_coordinate(e1,2))^2));
%             end
%         end
%     end
% end
% H = sparse(iH,jH,sH);
% Hs = sum(H,2);
%% INITIALIZE ITERATION
volfrac=0.4; penal=3; ft=2;
x = volfrac*ones(size(DESIGN_ZONE_ELEMENT_NODES,1),1);
xPhys = x(:);
% xPhys(x1>=3000&x4<=5000&z1<=min(z1)+dz)=0;
loop = 0;
change = 1;
m = 1;
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
c       = 10000*eeem;
d       = eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 1000;
kkttol  = 100;
load('gamma_structure')
load('U_tilder')
%
%%%% The iterations start:
kktnorm = kkttol+1000;
outit = 0;
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
F=[LOAD_MATRIX_def(:,1),Recovery_Matrix_def.',FAN_Center_Recovery_vector_def.'];
Res=K(free_dofs,free_dofs)\F(free_dofs,:);
U(free_dofs)=Res(:,1);
dK_dxU=zeros(length(U),length(xPhys(:)));
for xdim=1:length(xPhys(:))
    dK_dx=sparse(I(el_ij==xdim),J(el_ij==xdim),dK(el_ij==xdim),size(K_def,1),size(K_def,2)); dK_dx=dK_dx+triu(dK_dx.',1);
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
% if ft == 1
%     dperfo(:) = H*(x(:).*dperfo.')./Hs./max(1e-3,x(:));
% elseif ft == 2
%     dperfo = H*(dperfo.'./Hs);
%     dv(:) = H*(dv(:)./Hs);
% end
f0val=100*perfo;
fval=[sum(xPhys(:).*A2/2) - volfrac*sum(A2/2)];
df0dx=100*reshape(dperfo,[],1);
dfdx=[dv(:)'];
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
while change > 0.01
    outit   = outit+1;
    outeriter = outeriter+1;
    loop = loop + 1;
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    %     l1 = 0; l2 = 1e9; move = 0.2;
    %     while (l2-l1)/(l1+l2) > 1e-3
    %         lmid = 0.5*(l2+l1);
    %         xnew = max(0,max(x(:),min(1,min(x(:)+move,x(:).*sqrt(-dperfo.*(dperfo<0)./dv/lmid)))));
    %         if ft == 1
    %             xPhys = xnew;
    %         elseif ft == 2
    %             xPhys(:) = (H*xnew(:))./Hs;
    %         end
    %         if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
    %     end
    move = 0.25;
    lambda=dperfo;
    lambdav=-lambda/norm(lambda);
    mu=dv/norm(dv)*(sum(xPhys(:).*A2/2) - volfrac*sum(A2/2));
    muv=-mu/mean(A2);
    dir=lambdav(:)/norm(lambdav)-(muv/norm(muv))'*lambdav(:)/norm(lambdav)*(muv/norm(muv));
    xnew=max(0,min(1,x(:)+move*(dir(:))+muv));
%      if ft == 1
        xPhys = xnew;
%     elseif ft == 2
%         xPhys(:) = (H*xnew(:))./Hs;
%     end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    xval=x;
%     xval(x1>=3000&x4<=5000&z1<=min(z1)+dz)=0;
    
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
    F=[LOAD_MATRIX_def(:,1),Recovery_Matrix_def.',FAN_Center_Recovery_vector_def.'];
    Res=K(free_dofs,free_dofs)\F(free_dofs,:);
    U(free_dofs)=Res(:,1);
    dK_dxU=zeros(length(U),length(xPhys(:)));
    for xdim=1:length(xPhys(:))
        dK_dx=sparse(I(el_ij==xdim),J(el_ij==xdim),dK(el_ij==xdim),size(K_def,1),size(K_def,2)); dK_dx=dK_dx+triu(dK_dx.',1);
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
%     if ft == 1
%         dperfo(:) = H*(x(:).*dperfo.')./Hs./max(1e-3,x(:));
%     elseif ft == 2
%         dperfo = H*(dperfo.'./Hs);
%         dv(:) = H*(dv(:)./Hs);
%     end
%     xPhys(x1>=3000&x4<=5000&z1<=min(z1)+dz)=0;
%     dperfo(x1>=3000&x4<=5000&z1<=min(z1)+dz)=0;
    f0val=100*perfo;
    fval=[sum(xPhys(:).*A2/2) - volfrac*sum(A2/2)];
    df0dx=100*reshape(dperfo,[],1);
    dfdx=[dv(:)'];
    
    %% PLOT DENSITIES
    
    %     colormap(gray);h = figure(1); set(h,'Color',[1 1 1]);
    h = figure(1); set(h,'Color',[1 1 1]);
    hold on; h24_patch = patch('Vertices',DESIGN_ZONE_COORD(:,[2,4]),'Faces',DESIGN_ZONE_ELEMENT_NODES,'FaceVertexCData',(1-xval)*[1 1 1],'FaceColor','flat'); axis equal; axis off; drawnow;hold off
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',outeriter,perfo, ...
        sum(xPhys(:).*A2/2)/sum(A2/2),change);
    figure(2)
    hold on
    scatter(outeriter,f0val,'fill','k')
    figure(3)
    hold on
    scatter(outeriter,FAN_Ax_disp,'fill','b')
end
Final_coord=zeros(size(DESIGN_ZONE_COORD));
for k=1:size(DESIGN_ZONE_COORD,1)
    Final_coord(k,:)=DESIGN_ZONE_COORD(k,:)+[0,U(2*(k-1)+1),0,U(2*k)];
end
figure(4)
hold on; h24_patch = patch('Vertices',Final_coord(:,[2,4]),'Faces',DESIGN_ZONE_ELEMENT_NODES,'FaceVertexCData',(1-xval)*[1 1 1],'FaceColor','flat'); axis equal; axis off; drawnow;hold off
