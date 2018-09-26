%% ************* truss model code with MNA (curved) method ************* %%
clear ; close all ;
%% DATA
nelx = 80 ;
nely = 40 ;
penal = 3 ;
rmin = 1.2 ;
ft = 1 ;
C0 = 60 ;
T0 = 100 ;
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
P=4;
Sl=1;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
D=E0/(1-nu^2)*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
B=1/2*[-1 0 1 0 1 0 -1 0;0 -1 0 -1 0 1 0 1;-1 -1 -1 1 1 1 1 -1];
DB=D*B;
Cvm=[1 -0.5 0;-0.5 1 0;0 0 3];
Sel=DB'*Cvm*DB;
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
fixeddofs = [2*(fix(nelx/2)*(nely+1)+1:(nely+1):((nelx)*(nely+1)+1))-1;2*(fix(nelx/2)*(nely+1)+1:(nely+1):((nelx)*(nely+1)+1))];
fixeddofs=fixeddofs(:);

%% Generate the truss stiffness matrix
A=1;
I=100;
[Kcc,Kce,Kee,Fe,Fc,Recovery_matrixc,Recovery_matrixe]=truss_stiffness_no_condensation(nelx,nely,E0,A,I);
%% Generate the Projection matrix, load vector and stiffness matrix ready for assembly
coupling_nodes=nely+1:nely+1:(nely+1)*(nelx+1);
coupling_dofs=[coupling_nodes*2-1;coupling_nodes*2];
coupling_dofs=coupling_dofs(:);
Pr=sparse(coupling_dofs,1:2*(nelx+1),ones(1,2*(nelx+1)),2*(nelx+1)*(nely+1),2*(nelx+1));
Kt=[Pr*Kcc*Pr' Pr*Kce;Kce'*Pr' Kee];
Kt=(Kt+Kt')/2;
Ft=[Pr*Fc;Fe];
Rt=[Recovery_matrixc*Pr' zeros(size(Recovery_matrixc,1),size(Recovery_matrixe,2))
    zeros(size(Recovery_matrixe,1),size(Pr',2)) Recovery_matrixe];
U = zeros(2*(nely+1)*(nelx+1)+length(Fe),1);
F=U;
Lambda=[U,U];
F=F+Ft;
alldofs = 1:length(F);
freedofs = setdiff(alldofs,fixeddofs);
% original_force=-1;
Gamma=sparse(reshape( repmat(1:(nelx+1),2,1),[],1),reshape([1:(nelx+1);(1:(nelx+1))+(nelx+1)],[],1),reshape([ones(1,nelx+1);-ones(1,nelx+1)],[],1));
Rt=Gamma*Rt;
%% INITIALIZE ITERATION
% Starting configuration
nX = 12;             % Number of deformable elements along x
nY = 6;             % Number of deformable elements along y
lc = sqrt(nelx*nely/(0.25*nX*nY));
[xElts,yElts] = meshgrid(linspace(1/(nX+1)*nelx,nX/(nX+1)*nelx,nX),...
    linspace(1/(nY+1)*nely,nY/(nY+1)*nely,nY));
x = [reshape(xElts,1,numel(xElts));
    reshape(yElts,1,numel(yElts));
    zeros(1,numel(xElts));
    lc*ones(2,numel(xElts));0*ones(1,numel(xElts))];
x = reshape(x,numel(x),1);
n_tot = length(x);
n_var = 6 ;
Xmin    = repmat([0;0;-pi;0;0;-0.1],length(x)/6,1);
Xmax    = repmat([nelx;nely;pi;nelx;nely;0.1],length(x)/6,1);
x = (x-Xmin)./(Xmax-Xmin) ;

m = 3;
n = length(x);
eeen    = ones(n,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = x(:);
xold1   = xval;
xold2   = xval;
xmin    = zeron;
xmax    = eeen;
low     = xmin;
upp     = xmax;
C       = 1000*eeem;
D       = 0*eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 300;
kkttol  =0.01;
%
%%%% The iterations start:
kktnorm = kkttol+10;
% kktnorm = kkttol;
outit = 0;
[ud,Sd,~]=svds(Sel,8);
X_old=zeros(size(x));
X=Xmin+(Xmax-Xmin).*x;
%% START ITERATION
while kktnorm > kkttol && outit < maxoutit
    outit   = outit+1;
    outeriter = outeriter+1;
    %% density and its gradient calculation
    d = zeros(nelx*nely,1); % density
    dd = zeros(nelx*nely,n_tot); % gradient of density
    for in=1:n_tot/n_var  % loop on the components
        q = n_var*(in-1)+1;
        ln = sqrt(X(q+3)^2+X(q+4)^2);
        % FIND NEIGHBOURS
        limx = [max(min(floor(X(q)-ln),nelx),1) ; max(min(floor(X(q)+ln),nelx),1)];
        limy = [max(min(floor(X(q+1)-ln),nely),1) ; max(min(floor(X(q+1)+ln),nely),1)];
        [xN,yN] = meshgrid((limx(1):limx(2))-0.5,(limy(1):limy(2))-0.5);
        xN = reshape(xN,1,numel(xN));
        yN = reshape(yN,1,numel(yN));
        [ksi,eta] = changebasis3(xN,yN,X(q:q+n_var-1)) ;
        neigh = [ksi';eta'] ;
        nc = and(abs(neigh(1,:))<1,abs(neigh(2,:))<1);
        ind = reshape((ones(limy(2)-limy(1)+1,1)*(limx(1):limx(2))-1)*nely + ...
            (limy(1):limy(2))'*ones(1,limx(2)-limx(1)+1),(limx(2)-limx(1)+1)*(limy(2)-limy(1)+1),1);
        ind = ind(nc); % Global index of neighbours
        % COMPUTE DENSITIES
        if(~isempty(ind))
            % shape function and derivatives calculation
            v = neigh(:,nc) ;
            p1 = v(1,:)<0.5;
            p2 = v(2,:)<0.5;
            wksi = zeros(size(v,2),1);
            weta = zeros(size(v,2),1);
            dwdksi = zeros(size(v,2),1);
            dwdeta = zeros(size(v,2),1);
            wksi(p1) = 1-6*v(1,p1).^2+6*v(1,p1).^3;
            wksi(~p1) = 2-6*v(1,~p1)+6*v(1,~p1).^2-2*v(1,~p1).^3;
            weta(p2) = 1-6*v(2,p2).^2+6*v(2,p2).^3;
            weta(~p2) = 2-6*v(2,~p2)+6*v(2,~p2).^2-2*v(2,~p2).^3;
            f = wksi.*weta ;
            dwdksi(p1) = (-12*v(1,p1)+18*v(1,p1).^2);
            dwdksi(~p1) = (-6+12*v(1,~p1)-6*v(1,~p1).^2);
            dwdeta(p2) = (-12*v(2,p2)+18*v(2,p2).^2);
            dwdeta(~p2) = (-6+12*v(2,~p2)-6*v(2,~p2).^2);
            d(ind) = d(ind) + f;
            xN = xN(nc) ;
            yN = yN(nc) ;
            [dksidx,detadx] = grad_ksi_eta2(xN,yN,X(q:q+5),n_var) ;
            for j = 1:6
                dd(ind,q+j-1) = dwdksi.*dksidx(j,:)'.*wksi + dwdeta.*detadx(j,:)'.*weta ;
            end
        end
    end
    % ASYMPTOTIC DENSITY
    dd(d<1,:) = diag((1+12*d(d<1).^2-28*d(d<1).^3+15*d(d<1).^4))*dd(d<1,:);
    dd(d>=1,:) = 0;
    d = d+4*d.^3-7*d.^4+3*d.^5;
    d = sparse(min(d,1));
    dd = sparse(dd);
    %% FE-ANALYSIS
    
    sK = reshape(KE(:)*(Emin+d'.^penal*(E0-Emin)),64*nelx*nely,1);
    K0 = sparse(iK,jK,sK,size(Kt,1),size(Kt,2)); K=K0+Kt;K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
%      cndK=log(condest(K(freedofs,freedofs)))/log(10);
    
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    dv = ones(nely,nelx);
    %total compliance
    compliance=F'*U;
    ce_comp=reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    dcompliance = -penal*(E0-Emin)*dd'*(reshape(reshape(d,nely,nelx).^(penal-1).*ce_comp,nelx*nely,1));
    %relative squared displacement
%     rd1=(Gamma*Rt*U);

  mse = reshape(sqrt(sum(((U(edofMat)*ud(:,1:3))*Sd(1:3,1:3)).*(U(edofMat)*ud(:,1:3)),2)),nely,nelx); %microscopic Von Mises Stress
  rcv=reshape(d,nely,nelx).*(mse/Sl-1);
  G_ks=max(rcv(:))+1/P*log(mean(exp(P*(rcv(:)-max(rcv(:))))));
  rcv=d.*(mse(:)/Sl-1);
  sS=reshape(Sel(:)*(d'./mse(:)'/Sl.*exp(P*(rcv(:)'-max(rcv(:)))))/sum(exp(P*(rcv(:)-max(rcv(:))))),64*nelx*nely,1);
  S0=sparse(iK,jK,sS,length(F),length(F)); S0=(S0+S0')/2;
  dG_ksdu=S0*U;
%   Lambda2=U;
  
    rd=(Rt*U);
    c=((rd')*rd)/length(rd);
%     c1=((rd1')*rd1)/length(rd1);
  Lambda(freedofs,:)=K(freedofs,freedofs)\[(2*Rt(:,freedofs)'*rd),dG_ksdu(freedofs)];
  Lambda1=Lambda(:,1);
  Lambda2=Lambda(:,2);
%   Lambda2(freedofs) = K(freedofs,freedofs)\;
  dG_ksdx=dd'*((mse(:)/Sl-1).*exp(P*(rcv(:)-max(rcv(:))))/sum(exp(P*(rcv(:)-max(rcv(:))))));
  dG_du=reshape((sum((U(edofMat)*KE).*Lambda2(edofMat),2)),nely,nelx);
  dG_du =- penal*(E0-Emin)*dd'*(d.^(penal-1).*dG_du(:));
  dGKSl=dG_du(:)+dG_ksdx;
    ce=reshape(sum((Lambda1(edofMat)*KE).*U(edofMat),2),nely,nelx);
    dc = -1/length(rd)*penal*(E0-Emin)*dd'*(reshape(reshape(d,nely,nelx).^(penal-1).*ce,nelx*nely,1));
    %% Gather info for MMA
    f0val=mean(d)*100;
    fval=[(c-T0)/T0*100;G_ks*100;(compliance-C0)/C0*100];%;1*(-0.38-FAN_Ax_disp)/0.38
    fval=fval/1000;
    df0dx=mean(dd)'*100;
    dfdx=[dc(:)'/T0*100;dGKSl(:)'*100;dcompliance(:)'/C0*100];%;-1*dfandisp(:)'/0.38
    dfdx=dfdx/1000;
    %% MMA code optimization
    [x(:),ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,C,D);
    xold2 = xold1;
    xold1 = xval;
    xval  = x(:);
    change=norm(xval-xold1);
    %     xval(x1>=2500&x4<=4500&z1==min(z1))=0;
    
    %% %% The residual vector of the KKT conditions is calculated:
    [~,kktnorm,~] = ...
        kktcheck(m,n,x(:),ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,C,D);
        %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f kktnorm.:%7.3f ch.:%7.3f\n',outit,c, ...
        full(f0val),kktnorm,change);
        %% PLOT DENSITIES
      figure(1)
%       colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    h = figure(1); set(h,'Color',[1 1 1]);
    [Yy,Xx]=find(nodenrs);
    Yy=nely+1-Yy;
    Xx=Xx-1;
    hold on;  patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-d)*[1 1 1],'FaceColor','flat','EdgeColor','none'); axis equal; axis off; drawnow;hold off
    title(['Design zone at iteration ',num2str(outit)])
    colormap(gray);
    % PLOT MEMBERS
    X=Xmin+(Xmax-Xmin).*x;
    figure(2)
    fact = 0.25 ; % entre 0.2 et 0.3 
    [gridx,gridy] = meshgrid(0:nelx,0:nely) ;
    Phi=cell(n_tot/n_var,1);
    for i=1:n_tot/n_var
        q = 6*(i-1)+1;
        Phi{i} = Phi_c(X(q:q+5),gridx,gridy,6,fact) ;
    end
    Phi_max = Phi{1} ;
    for j = 2:n_tot/n_var
        Phi_max = max(Phi_max,Phi{j}) ;
    end
    contourf(gridx,gridy,Phi_max,[0,0]) ;
    axis equal; axis([0 nelx 0 nely]);  drawnow;
end