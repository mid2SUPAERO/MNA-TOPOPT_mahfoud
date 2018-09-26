function [Kcc,Kce,Kee,Fe,Fc,Recovery_matrixc,Recovery_matrixe,Fte,Ftc,X,Y,DX,DY,ELEMENT]=truss_stiffness_no_condensation(nelx,nely,E,A,I)
L=1;
Kt=E*A/L;
Kf=E*I/L^3;
Keo=[Kt 0 0 -Kt 0 0
    0 12*Kf 6*L*Kf 0 -12*Kf 6*L*Kf
    0 6*L*Kf 4*L^2*Kf 0 -6*L*Kf 2*L^2*Kf
    -Kt 0 0 Kt 0 0
    0 -12*Kf -6*L*Kf 0 12*Kf -6*L*Kf
    0 6*L*Kf 2*L^2*Kf 0 -6*L*Kf 4*L^2*Kf];
Pov=[0 -1 0 0 0 0
    1 0 0 0 0 0
    0 0 1 0 0 0
    0 0 0 0 -1 0
    0 0 0 1 0 0
    0 0 0 0 0 1];
Kev=Pov'*Keo*Pov;
% node_list=1:(2*(nelx+1)+3*(nely-1));
node_dofs_map=reshape((1:3*(2*(nelx+1)+3*(nely-1)))',3,[])';
ELEMENT(1:nelx,:)=[(1:nelx)',(1:nelx)'+1];
ELEMENT(nelx+(1:nelx),:)=[nelx+1+(1:nelx)',(1:nelx)'+nelx+2];
ELEMENT(2*nelx+(1:(nely)),:)=[1,2*nelx+2+1
                              ((2*nelx+2+1):(2*nelx+2+nely-2))',((2*nelx+2+1):(2*nelx+2+nely-2))'+1
                                                    (2*nelx+2+nely-1), nelx+2];
ELEMENT(2*nelx+nely+(1:(nely)),:)=[fix(nelx/2),2*nelx+2+nely
                                                    ((2*nelx+2+nely):(2*nelx+2+2*nely-3))',((2*nelx+2+nely):(2*nelx+2+2*nely-3))'+1
                                                    (2*nelx+2+2*nely-2), nelx+1+fix(nelx/2)];
ELEMENT(2*nelx+2*nely+(1:(nely)),:)=[nelx+1,2*nelx+2+2*nely-2+1
                                                    ((2*nelx+2+2*nely-2+1):(2*nelx+2+3*(nely-1)-1))',((2*nelx+2+2*nely-2+1):(2*nelx+2+3*(nely-1)-1))'+1
                                                    (2*nelx+2+3*(nely-1)), 2*nelx+2];
ELEMENT_DOF=[node_dofs_map(ELEMENT(:,1),:),node_dofs_map(ELEMENT(:,2),:)];
[io,jo,ko]=find(Keo);
[iv,jv,kv]=find(Kev);
Io=reshape(ELEMENT_DOF(1:2*nelx,io)',[],1);
Jo=reshape(ELEMENT_DOF(1:2*nelx,jo)',[],1);
% Ko=repmat(ko,2*nelx,1);
Ko=[repmat(ko,nelx,1)/4;4*repmat(ko,nelx,1)];
%% thermal equivalent loads
Bo=[-1/L 0 0 1/L 0 0];
Bv=Bo*Pov;
alphat=12e-6;
DT=200;
fo=E*A*alphat*DT*L*Bo;
[ifo,~,kfo]=find(fo(:));
Ifo=reshape(ELEMENT_DOF(1:2*nelx,ifo)',[],1);
Kfo=[repmat(kfo(:),nelx,1)/4;4*repmat(kfo(:),nelx,1)];
fv=E*A*alphat*DT*L*Bv;
[ifv,~,kfv]=find(fv(:));
Ifv=reshape(ELEMENT_DOF(2*nelx+(1:3*nely),ifv)',[],1);
Kfv=4*repmat(kfv(:),3*nely,1);
Ft=sparse([Ifo;Ifv],ones(size([Ifo;Ifv])),[Kfo;Kfv],max(node_dofs_map(:)),1);
Iv=reshape(ELEMENT_DOF(2*nelx+(1:3*nely),iv)',[],1);
Jv=reshape(ELEMENT_DOF(2*nelx+(1:3*nely),jv)',[],1);
Kv=4*repmat(kv,3*nely,1);
K=sparse([Io;Iv],[Jo;Jv],[Ko;Kv],max(node_dofs_map(:)),max(node_dofs_map(:)));
K=(K+K')/2;
F=sparse(3*((nelx+1))+1,1,-1,max(node_dofs_map(:)),1);
coupling_dofs=[(1:3:(3*(nelx+1)-2));2:3:(3*(nelx+1)-1)];
coupling_dofs=coupling_dofs(:);
reduced_dofs=setdiff(node_dofs_map(:),coupling_dofs);
% total compliance
% observation_dofs=[3*((nelx+1))+1];
%squared relative displacement
observation_dofs=[2:3:3*((nelx+1));(2+3*((nelx+1))):3:3*(2*(nelx+1))];
observation_dofs=observation_dofs(:)';
u0=zeros(size(K,1),1);
u0(reduced_dofs)=K(reduced_dofs,reduced_dofs)\F(reduced_dofs);
Kcc=K(coupling_dofs,coupling_dofs);
Kce=K(coupling_dofs,reduced_dofs);
Kee=K(reduced_dofs,reduced_dofs);
Fc=F(coupling_dofs);
Fe=F(reduced_dofs);
Ftc=Ft(coupling_dofs);
Fte=Ft(reduced_dofs);
% B=(Kce/Kee);
% Kcctilde=Kcc-B*(Kce');
% Ftilde=Fc-B*Fe;
% total compliance
%Recovery_matrix=-sparse(1:length(observation_dofs),find(Fe),ones(size(observation_dofs)),length(observation_dofs),length(reduced_dofs))*B';
%squared relative displacement
selector=zeros(size(u0));
selector(observation_dofs)=1;
selector1=sparse(selector(reduced_dofs));
selector2=sparse(selector(coupling_dofs));
% Recovery_matrix=[-sparse(1:length(find(selector1)),find(selector1),ones(size(find(selector1))),length(find(selector1)),length(reduced_dofs))*B';...
%     sparse(1:length(find(selector2)),find(selector2),ones(size(find(selector2))),length(find(selector2)),length(coupling_dofs))];
Recovery_matrixc=sparse(1:length(find(selector2)),find(selector2),ones(size(find(selector2))),length(find(selector2)),length(coupling_dofs));
Recovery_matrixe=sparse(1:length(find(selector1)),find(selector1),ones(size(find(selector1))),length(find(selector1)),length(reduced_dofs));
Xo=[(0:nelx)';(0:nelx)'];
Yo=[zeros(nelx+1,1);repmat(-nely,nelx+1,1)];
Xv=[zeros(nely-1,1);repmat(fix(nelx/2)-1,nely-1,1);repmat(nelx,nely-1,1)];
Yv=repmat((-1:-1:-nely+1)',3,1);
X=[Xo;Xv];
Y=[Yo;Yv];
plot(X(ELEMENT)',Y(ELEMENT)','k-o','MarkerFaceColor','k')
axis equal
U=reshape(u0,3,[])';
U=U(:,1:2);
DX=U(:,1);
DY=U(:,2);
hold on
plot(X(ELEMENT)'+DX(ELEMENT)',Y(ELEMENT)'+DY(ELEMENT)','-.r')
% u0=u0(observation_dofs);
