%% ***************** MNA code with curved components ******************* %%
clear ; close all ; tic
%% data
nelx = 8*10;          % Number of finite elements along x
nely = 8*10;          % Number of finite elements along y
volfrac = 0.4*0.75 ;     % Targeted Volume fraction
penal0 = 3 ;
% Starting configuration
nX = 12 ;             % Number of deformable elements along x
nY = 12 ;             % Number of deformable elements along y
tolchange=1e-3 ;
%Creation of the initial x0 (regular grid of moving node)
d = sqrt(volfrac*nelx*nely/(0.5*nX*nY));
[xElts,yElts] = meshgrid(linspace(1/(nX+1)*nelx,nX/(nX+1)*nelx,nX),...
    linspace(1/(nY+1)*nely,nY/(nY+1)*nely,nY));
% x = [reshape(xElts,1,numel(xElts));
%     reshape(yElts,1,numel(yElts));
%     zeros(1,numel(xElts));
%     d*ones(2,numel(xElts));0.01*ones(1,numel(xElts))];
% L-shape
rm = and(xElts > 0.5*nelx,yElts > 0.5*nely);
xElts(rm) = nan;
yElts(rm) = nan;
x = [reshape(xElts(~isnan(xElts)),1,numel(xElts(~isnan(xElts))));
    reshape(yElts(~isnan(yElts)),1,numel(yElts(~isnan(yElts))));
    zeros(1,numel(xElts(~isnan(xElts))));
    d*ones(2,numel(xElts(~isnan(xElts))));0.01*ones(1,numel(xElts(~isnan(xElts))))];
x = reshape(x,numel(x),1);
n_var = 6 ;
n_tot = length(x);
% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
P = 20; Sl = 0.5 ;
rhoMax = 1.05;
mMax = nelx*nely*volfrac;
penal = 3;
dp = (penal0-penal)/5; % initial penalization
f1 = figure('units','normalized','position',[0.1,0.25,0.3,0.5]);
f2 = figure('units','normalized','position',[0.6,0.25,0.3,0.5]);
f3 = figure() ;
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
% DEFINE LOADS AND SUPPORTS (Cantilever,L-shpae or MBB beam)
% F = sparse(2*(nely+1)*nelx+nely,1,-1,2*(nely+1)*(nelx+1),1); % Cantilever
% F = sparse(2*(nely+1),1,-1,2*(nely+1)*(nelx+1),1); % MBB beam
F = sparse(2*(nely+1)*nelx+nely,1,-1,2*(nely+1)*(nelx+1),1); % L-shape
U = zeros(2*(nely+1)*(nelx+1),1);
% fixeddofs = 1:2*(nely+1); % Cantilever
% fixeddofs = union(1:2:2*(nely+1),2*nelx*(nely+1)+2); % MBB beam
fixeddofs = kron(2*nely:2*(nely+1):(nelx+2)*(nely+1),[1,1])...
    + repmat([1,2],1,nelx/2+1); % L-shape
% emptyelts = [] ; % Cantilever and MBB
emptyelts = repmat(0.5*nely*(nelx+1)+1:1:nely*(0.5*nelx+1),1,0.5*nelx)...
    +kron(0:nely:(0.5*nelx-1)*nely,ones(1,0.5*nely)); % L-shape
fullelts = [];
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
% normalization of variables
Xmin    = repmat([0;0;-pi;0;0;-0.05],n_tot/n_var,1);
Xmax    = repmat([nelx;nely;pi;nelx;nely;0.05],n_tot/n_var,1);
x = (x-Xmin)./(Xmax-Xmin) ;
%% INITIALIZE ITERATION
change = 1;
mm = 2 ;
eeen    = ones(n_tot,1);
eeem    = ones(mm,1);
zerom   = zeros(mm,1);
xval    = x(:);
xold1   = xval;
xold2   = xval;
xmin    = 0*eeen;
xmax    = eeen;
low     = xmin;
upp     = xmax;
C       = 1000*eeem;
d_       = 0*eeem;
a0      = 1;
aa       = zerom;
outeriter = 0;
maxoutit  = 600;
kkttol  =tolchange;
kktnorm = kkttol+10;
outit = 0;
[ud,Sd,~]=svds(Sel,8);
X_old=zeros(size(x));
X=Xmin+(Xmax-Xmin).*x;
%% START ITERATION
while (kktnorm > kkttol && outit < maxoutit)&&change>tolchange
    outit   = outit+1;
    outeriter = outeriter+1;
    delta = 1e-3 ;
    dd = zeros(nelx*nely,n_tot) ;
    %     [~,~,~,~,d] = FEA_comp_vol(X) ;
    d = zeros(nelx*nely,1); % density
    for in=1:n_tot/n_var  % loop on the components
        n = n_var*(in-1)+1;
        ln = sqrt(X(n+3)^2+X(n+4)^2);
        % FIND NEIGHBOURS
        limx = [max(min(floor(X(n)-ln),nelx),1) ; max(min(floor(X(n)+ln),nelx),1)];
        limy = [max(min(floor(X(n+1)-ln),nely),1) ; max(min(floor(X(n+1)+ln),nely),1)];
        [xN,yN] = meshgrid((limx(1):limx(2))-0.5,(limy(1):limy(2))-0.5);
        xN = reshape(xN,1,numel(xN));
        yN = reshape(yN,1,numel(yN));
        [ksi,eta] = changebasis3(xN,yN,X(n:n+n_var-1)) ;
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
            wksi(p1) = 1-6*v(1,p1).^2+6*v(1,p1).^3;
            wksi(~p1) = 2-6*v(1,~p1)+6*v(1,~p1).^2-2*v(1,~p1).^3;
            weta(p2) = 1-6*v(2,p2).^2+6*v(2,p2).^3;
            weta(~p2) = 2-6*v(2,~p2)+6*v(2,~p2).^2-2*v(2,~p2).^3;
            f = wksi.*weta ;
            d(ind) = d(ind) + f;
        end
    end
    for k = 1:n_tot
        pert = zeros(n_tot,1) ;
        pert(k) = delta ;
        x1 = x-pert ;
        x2 = x+pert ;
        X1 = Xmin+(Xmax-Xmin).*x1 ;
        X2 = Xmin+(Xmax-Xmin).*x2 ;
        %         [~,~,~,~,d1] = FEA_comp_vol(X1) ;
        %         [~,~,~,~,d2] = FEA_comp_vol(X2) ;
        d1 = zeros(nelx*nely,1); % density
        for in=1:n_tot/n_var  % loop on the components
            n = n_var*(in-1)+1;
            ln = sqrt(X1(n+3)^2+X1(n+4)^2);
            % FIND NEIGHBOURS
            limx = [max(min(floor(X1(n)-ln),nelx),1) ; max(min(floor(X1(n)+ln),nelx),1)];
            limy = [max(min(floor(X1(n+1)-ln),nely),1) ; max(min(floor(X1(n+1)+ln),nely),1)];
            [xN,yN] = meshgrid((limx(1):limx(2))-0.5,(limy(1):limy(2))-0.5);
            xN = reshape(xN,1,numel(xN));
            yN = reshape(yN,1,numel(yN));
            [ksi,eta] = changebasis3(xN,yN,X1(n:n+n_var-1)) ;
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
                wksi(p1) = 1-6*v(1,p1).^2+6*v(1,p1).^3;
                wksi(~p1) = 2-6*v(1,~p1)+6*v(1,~p1).^2-2*v(1,~p1).^3;
                weta(p2) = 1-6*v(2,p2).^2+6*v(2,p2).^3;
                weta(~p2) = 2-6*v(2,~p2)+6*v(2,~p2).^2-2*v(2,~p2).^3;
                f = wksi.*weta ;
                d1(ind) = d1(ind) + f;
            end
        end
        d2 = zeros(nelx*nely,1); % density
        for in=1:n_tot/n_var  % loop on the components
            n = n_var*(in-1)+1;
            ln = sqrt(X2(n+3)^2+X2(n+4)^2);
            % FIND NEIGHBOURS
            limx = [max(min(floor(X2(n)-ln),nelx),1) ; max(min(floor(X2(n)+ln),nelx),1)];
            limy = [max(min(floor(X2(n+1)-ln),nely),1) ; max(min(floor(X2(n+1)+ln),nely),1)];
            [xN,yN] = meshgrid((limx(1):limx(2))-0.5,(limy(1):limy(2))-0.5);
            xN = reshape(xN,1,numel(xN));
            yN = reshape(yN,1,numel(yN));
            [ksi,eta] = changebasis3(xN,yN,X2(n:n+n_var-1)) ;
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
                wksi(p1) = 1-6*v(1,p1).^2+6*v(1,p1).^3;
                wksi(~p1) = 2-6*v(1,~p1)+6*v(1,~p1).^2-2*v(1,~p1).^3;
                weta(p2) = 1-6*v(2,p2).^2+6*v(2,p2).^3;
                weta(~p2) = 2-6*v(2,~p2)+6*v(2,~p2).^2-2*v(2,~p2).^3;
                f = wksi.*weta ;
                d2(ind) = d2(ind) + f;
                dd(:,k) = (d2-d1)/(X2(k)-X1(k)) ;
            end
        end
    end
    % ASYMPTOTIC DENSITY
    dd(d<1,:) = diag((1+12*d(d<1).^2-28*d(d<1).^3+15*d(d<1).^4))*dd(d<1,:);
    dd(d>=1,:) = 0;
    d = d+4*d.^3-7*d.^4+3*d.^5;
    d(emptyelts) = 0; dd(emptyelts,:) = 0;
    d(fullelts) = 1; dd(fullelts,:) = 0;
    d = sparse(min(d,1));
    dd = sparse(dd);
    % FEM analysis
    sK = reshape(KE(:)*(Emin+d'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    % compliance and it's gradient calculation
    ce = sum((U(edofMat)*KE).*U(edofMat),2);
    c = sum(sum((Emin+d.^penal*(E0-Emin)).*ce));
    dc = -penal*(E0-Emin)*dd'*(d.^(penal-1).*ce);
    % volume constraint
    cst=(sum(d(:))-mMax)/mMax;
    dcst=sum(dd)/mMax;
    % stress constraint
    mse = reshape(sqrt(sum(((U(edofMat)*ud(:,1:3))*Sd(1:3,1:3)).*(U(edofMat)*ud(:,1:3)),2)),nely,nelx); %microscopic Von Mises Stress
    rcv=reshape(d,nely,nelx).*(mse/Sl-1);
    G_ks=max(rcv(:))+1/P*log(mean(exp(P*(rcv(:)-max(rcv(:))))));
    rcv=d.*(mse(:)/Sl-1);
    sS=reshape(Sel(:)*(d'./mse(:)'/Sl.*exp(P*(rcv(:)'-max(rcv(:)))))/sum(exp(P*(rcv(:)-max(rcv(:))))),64*nelx*nely,1);
    S0=sparse(iK,jK,sS,length(F),length(F)); S0=(S0+S0')/2;
    dG_ksdu=S0*U;
    Lambda=U;
    Lambda(freedofs) = K(freedofs,freedofs)\dG_ksdu(freedofs);
    dG_ksdx=dd'*((mse(:)/Sl-1).*exp(P*(rcv(:)-max(rcv(:))))/sum(exp(P*(rcv(:)-max(rcv(:))))));
    dG_du=reshape((sum((U(edofMat)*KE).*Lambda(edofMat),2)),nely,nelx);
    dG_du =- penal*(E0-Emin)*dd'*(d.^(penal-1).*dG_du(:));
    dGKSl=dG_du(:)+dG_ksdx;
    dc = dc.*(Xmax-Xmin) ; dcst = dcst.*((Xmax-Xmin)') ; dGKSl = dGKSl.*(Xmax-Xmin) ;
    dc=dc.*(Xmax-Xmin);dcst=dcst.*((Xmax-Xmin)');
    change=max(abs(X_old-X));
    X_old=X;
    f0val=c/10000;
    fval=[cst;G_ks*10] ;
    df0dx=dc(:)/10000;
    dfdx_=[dcst(:)';dGKSl(:)'*10] ;
    % MMA code optimization
    [x(:),ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(mm,n_tot,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx_,low,upp,a0,aa,C,d_);
    xold2 = xold1;
    xold1 = xval;
    xval  = x(:);
    % The residual vector of the KKT conditions is calculated:
    [~,kktnorm,~] = ...
        kktcheck(mm,n_tot,x(:),ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx_,a0,aa,C,d_);
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f Stress.:%7.3f kktnorm.:%7.3f ch. %7.3f\n',outit,c, ...
        cst,full(G_ks),kktnorm,change);
    % PLOT MEMBERS
    X=Xmin+(Xmax-Xmin).*x;
    figure(f1)
    fact = 0.4 ; % entre 0.2 et 0.3
    [gridx,gridy] = meshgrid(0:nelx,0:nely) ;
    Phi=cell(n_tot/n_var,1);
    for i=1:n_tot/n_var
        n = 6*(i-1)+1;
        Phi{i} = Phi_c(X(n:n+5),gridx,gridy,6,fact) ;
    end
    Phi_max = Phi{1} ;
    for j = 2:n_tot/n_var
        Phi_max = max(Phi_max,Phi{j}) ;
    end
    contourf(gridx,gridy,Phi_max,[0,0]) ;
    axis equal; axis([0 nelx 0 nely]);  drawnow; title('component plot') ;
    figure(f2)
    colormap(flipud(gray)); imagesc(flipud(reshape(d,nely,nelx))); caxis([0 1]); axis equal; axis off; drawnow;
    title('density field') ;
    if outit>20
        figure(f3)
        subplot(3,1,1)
        scatter(outit,c,4) ; xlabel('iteration') ; ylabel('compliance');
        title('compliance history') ;
        hold on
        subplot(3,1,2)
        scatter(outit,kktnorm,4) ; xlabel('iteration') ; ylabel('KKT norm');
        title('KKT norm history') ;
        hold on
        subplot(3,1,3)
        scatter(outit,cst,4) ; xlabel('iteration') ; ylabel('mass constraint (in %)');
        title('mass constraint history') ;
        hold on
    end
end
%% stress distribution
figure ;
colormap(hot) ; imagesc(flipud(reshape(rcv,nely,nelx))) ; 
caxis([min(rcv(:)) max(rcv(:))]) ; axis equal ; drawnow ; colorbar ;
title('relaxed constraint violation') ;
toc