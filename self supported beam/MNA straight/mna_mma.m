%% ***************** MNA code mma with old variables ******************* %%
clear ; close all ; tic
%% data
nelx = 16*10;          % Number of finite elements along x
nely = 8*10;          % Number of finite elements along y
volfrac = 0.5;     % Targeted Volume fraction
penal0 = 3 ;
% Starting configuration
nX = 8;             % Number of deformable elements along x
nY = 4;             % Number of deformable elements along y
tolchange=1e-3;
%Creation of the initial x0 (regular grid of moving node)
d = sqrt(volfrac*nelx*nely/(0.5*nX*nY));
[xElts,yElts] = meshgrid(linspace(1/(nX+1)*nelx,nX/(nX+1)*nelx,nX),...
    linspace(1/(nY+1)*nely,nY/(nY+1)*nely,nY));
x = [reshape(xElts,1,numel(xElts));
    reshape(yElts,1,numel(yElts));
    zeros(1,numel(xElts));
    d*ones(2,numel(xElts))];
x = reshape(x,numel(x),1);
comp_opt = 100 ; x_opt = x ;
% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
rhoMax = 1.05;
mMax = nelx*nely*volfrac ;
penal = 3 ;
dp = (penal0-penal)/5; % initial penalization
f1 = figure('units','normalized','position',[0.1,0.25,0.3,0.5]);
f2 = figure('units','normalized','position',[0.6,0.25,0.3,0.5]);
f3 = figure() ;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% define boundary conditions (self suppoted beam)
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = [2 2*(nely+1)*nelx+2] ; 
emptyelts = [] ;
fullelts = [];
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
% adimensional variable
Xmin    = repmat([0;0;-pi;0;0],length(x)/5,1);
Xmax    = repmat([nelx;nely;pi;nelx;nely],length(x)/5,1);
x = (x-Xmin)./(Xmax-Xmin) ;
%% saving pictures
% parent_dir_name = 'mesh80x40_old_var_v0.5\';
% mkdir(parent_dir_name);
%% INITIALIZE ITERATION
change = 1;
mm = 1;
n_ = length(x(:));
eeen    = ones(n_,1);
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
maxoutit  = 300;
kkttol  =tolchange;
kktnorm = kkttol+10;
outit = 0;
X_old=zeros(size(x));
X=Xmin+(Xmax-Xmin).*x;
%% START ITERATION
while (kktnorm > kkttol && outit < maxoutit)&&change>tolchange
    outit   = outit+1;
    outeriter = outeriter+1;
    % density and it's gradient calculation
    d = zeros(nelx*nely,1); % density
    dd = zeros(nelx*nely,length(X)); % gradient of density
    for in=1:length(X)/5  % loop on the components
        n = 5*(in-1)+1;
        ln = sqrt(X(n+3)^2+X(n+4)^2);
        sn = sin(X(n+2));
        cn = cos(X(n+2));
        % FIND NEIGHBOURS
        limx = [max(min(floor(X(n)-ln),nelx),1) ; max(min(floor(X(n)+ln),nelx),1)];
        limy = [max(min(floor(X(n+1)-ln),nely),1) ; max(min(floor(X(n+1)+ln),nely),1)];
        [xN,yN] = meshgrid((limx(1):limx(2))-0.5,(limy(1):limy(2))-0.5);
        xN = reshape(xN,1,numel(xN));
        yN = reshape(yN,1,numel(yN));
        neigh = [ ((xN-X(n))*cn + (yN-X(n+1))*sn)/X(n+3);
            (-(xN-X(n))*sn + (yN-X(n+1))*cn)/X(n+4)];
        nc = and(abs(neigh(1,:))<1,abs(neigh(2,:))<1);
        ind = reshape((ones(limy(2)-limy(1)+1,1)*(limx(1):limx(2))-1)*nely + ...
            (limy(1):limy(2))'*ones(1,limx(2)-limx(1)+1),(limx(2)-limx(1)+1)*(limy(2)-limy(1)+1),1);
        ind = ind(nc); % Global index of neighbours
        % COMPUTE DENSITIES
        if(~isempty(ind))
            % shape function and derivatives calculation
            v = neigh(:,nc) ;
            sx = sign(v(1,:)');
            sy = sign(v(2,:)');
            v = abs(v);
            
            p1 = v(1,:)<0.5;
            p2 = v(2,:)<0.5;
            
            wx = zeros(size(v,2),1);
            wy = zeros(size(v,2),1);
            dfdx = zeros(size(v,2),1);
            dfdy = zeros(size(v,2),1);
            
            wx(p1) = 1-6*v(1,p1).^2+6*v(1,p1).^3;
            wx(~p1) = 2-6*v(1,~p1)+6*v(1,~p1).^2-2*v(1,~p1).^3;
            wy(p2) = 1-6*v(2,p2).^2+6*v(2,p2).^3;
            wy(~p2) = 2-6*v(2,~p2)+6*v(2,~p2).^2-2*v(2,~p2).^3;
            
            f = abs(wx.*wy);
            
            dfdx(p1) = (-12*v(1,p1)+18*v(1,p1).^2);
            dfdx(~p1) = (-6+12*v(1,~p1)-6*v(1,~p1).^2);
            dfdx = dfdx.*sx.*wy;
            dfdy(p2) = (-12*v(2,p2)+18*v(2,p2).^2);
            dfdy(~p2) = (-6+12*v(2,~p2)-6*v(2,~p2).^2);
            dfdy = dfdy.*sy.*wx;
            
            dfdx = dfdx/X(n+3); dfdy = dfdy/X(n+4);
            d(ind) = d(ind) + f;
            dd(ind,n) = - cn*dfdx + sn*dfdy;
            dd(ind,n+1) = - sn*dfdx - cn*dfdy;
            dd(ind,n+2) = (-(xN(nc)-X(n))*sn + (yN(nc)-X(n+1))*cn)'.*dfdx - ...
                ((xN(nc)-X(n))*cn + (yN(nc)-X(n+1))*sn)'.*dfdy;
            dd(ind,n+3) = - ((xN(nc)-X(n))*cn + (yN(nc)-X(n+1))*sn)'.*dfdx/X(n+3);
            dd(ind,n+4) = ((xN(nc)-X(n))*sn - (yN(nc)-X(n+1))*cn)'.*dfdy/X(n+4);
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
    % FEM analysis & loads update
    F = sparse((nely+1)*nelx+nely,1,-sum(d(:))/mMax,2*(nely+1)*(nelx+1),1);
    sK = reshape(KE(:)*(Emin+d'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    % compliance and it's gradient calculation
    dg = -1/mMax*dd'*ones(nelx*nely,1) ;
    dF = zeros(2*(nelx+1)*(nely+1),n_) ;
    dF((nely+1)*nelx+nely,:) = dg' ; dF = sparse(dF) ;
    ce = sum((U(edofMat)*KE).*U(edofMat),2) ;
    c = sum(sum((Emin+d.^penal*(E0-Emin)).*ce)) ;
    dc = -penal*(E0-Emin)*dd'*(d.^(penal-1).*ce) + 2*dF'*U ;
    cst=-(sum(d(:))-mMax)/mMax*100;
    dcst=-sum(dd)/mMax*100;
    dc=dc.*(Xmax-Xmin);dcst=dcst.*((Xmax-Xmin)');
    change=max(abs(X_old-X));
    X_old=X;
    f0val=c/10000;
    fval=cst;
    df0dx=dc(:)/10000;
    dfdx_=dcst(:)';
    % MMA code optimization
    [x(:),ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(mm,n_,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx_,low,upp,a0,aa,C,d_);
    xold2 = xold1;
    xold1 = xval;
    xval  = x(:);
    fprintf(' It.:%5i Obj.:%11.4f Const.:%7.3f kktnorm.:%7.3f ch. %7.3f\n',outit,c, ...
        cst,kktnorm,change);
    % The residual vector of the KKT conditions is calculated:
    [~,kktnorm,~] = ...
        kktcheck(mm,n_,x(:),ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx_,a0,aa,C,d_);
    % PLOT MEMBERS
    X=Xmin+(Xmax-Xmin).*x;
    figure(f1)
    %     plotConfig(X);
    fact = 0.45;%.37
    clf
    hold on
    for i=1:length(X)/5
        n = 5*(i-1)+1;
        sn = sin(X(n+2));
        cn = cos(X(n+2));
        fill([X(n)+fact*X(n+3)*cn-fact*X(n+4)*sn,...
            X(n)-fact*X(n+3)*cn-fact*X(n+4)*sn,...
            X(n)-fact*X(n+3)*cn+fact*X(n+4)*sn,...
            X(n)+fact*X(n+3)*cn+fact*X(n+4)*sn,...
            X(n)+fact*X(n+3)*cn-fact*X(n+4)*sn],...
            [X(n+1)+fact*X(n+3)*sn+fact*X(n+4)*cn,...
            X(n+1)-fact*X(n+3)*sn+fact*X(n+4)*cn,...
            X(n+1)-fact*X(n+3)*sn-fact*X(n+4)*cn,...
            X(n+1)+fact*X(n+3)*sn-fact*X(n+4)*cn,...
            X(n+1)+fact*X(n+3)*sn+fact*X(n+4)*cn],[0 0.2 0.4])
    end
    hold off
    axis equal; axis([0 nelx 0 nely]); drawnow ; title('component plot') ;
    % save pictures
%     FileName=[parent_dir_name,'\Fig1_',int2str(outit),'.png'];
%     saveas(1,FileName);
    figure(f2)
    colormap(flipud(gray)); imagesc(flipud(reshape(d,nely,nelx))); caxis([0 1]); axis equal; axis off; drawnow;
    title('density field') ;
    % mise à jour de la solution optimale
    if cst<=0 && outit>20 && c<comp_opt 
        x_opt=x; comp_opt=c ;
    end
    % convergence history
    if outit>15
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
toc