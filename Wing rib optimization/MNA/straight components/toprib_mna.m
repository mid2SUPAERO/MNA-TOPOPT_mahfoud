%% ************ Wing rib topology optimization with MNA **************** %%
% this code presents a topology optimization of a wing rib under MNA
% framework with straight components
clear ; close all ;
%% data
nelx=160;
nely=40;
volfrac=0.45;
penal=3;
tol = 0.01 ;
E0 = 1;
Emin = 1e-9;
nu = 0.3;
mMax = nelx*nely*volfrac ;
d = (3*nelx)/8;
R = 6*nely;                                   % radii
Rint1 = nely/4;
Rint2 = nely/4;
Rint3 = nely/6;
a = floor((sqrt((R)^2-(d)^2))-R+nely);        % edges
b = floor((sqrt((R)^2-(nelx-d)^2))-R+nely);
f = ((nelx/4)+Rint1+(((nelx/4)-Rint1-Rint2)/2)-(nelx/80)+1);  %white square dimensions
g = ((nelx/4)+Rint1+(((nelx/4)-Rint1-Rint2)/2)+(nelx/80));
%% initial design
nX = 16 ; % number of components along x axis
nY = 4 ; % number of components along y axis
l_ini = sqrt(volfrac*nelx*nely/(0.5*nX*nY)) ; % initial length of component
[xElts,yElts] = meshgrid(linspace(1/(nX+1)*nelx,nX/(nX+1)*nelx,nX),...
    linspace(1/(nY+1)*nely,nY/(nY+1)*nely,nY)) ;
xElts = reshape(xElts,1,numel(xElts)) ;
yElts = reshape(yElts,1,numel(yElts)) ;
% % removing components
% for i = 1:length(xElts)
%     % white external round
%     if sqrt((yElts(i)-R)^2+(xElts(i)-d)^2) > R
%         xElts(i) = nan ; yElts(i) = nan ;
%     end
% end
% for i = 1:length(xElts)
%     % white internal round 1
%     if sqrt((yElts(i)-((3*nely)/5))^2+(xElts(i)-(nelx/4))^2) <= Rint1
%         xElts(i) = nan ; yElts(i) = nan ;
%     end
% end
% for i = 1:length(xElts)
%     % white internal round 2
%     if sqrt((yElts(i)-((3*nely)/5))^2+(xElts(i)-(nelx/2))^2) <= Rint2
%         xElts(i) = nan ; yElts(i) = nan ;
%     end
% end
% for i = 1:length(xElts)
%     % white internal round 3
%     if sqrt((yElts(i)-((3*nely)/5))^2+(xElts(i)-((3*nelx)/4))^2) <= Rint3
%         xElts(i) = nan ; yElts(i) = nan ;
%     end
% end
% for i = 1:length(xElts)
%     % white square
%     if f < xElts(i) < g && ((nely/10)+1) < yElts(i) < (9*nely)/10
%         xElts(i) = nan ; yElts(i) = nan ;
%     end
% end
xElts = xElts(~isnan(xElts)) ;
yElts = yElts(~isnan(yElts)) ;
n_c = length(xElts) ;
x = [xElts;yElts;zeros(1,n_c);l_ini*ones(2,n_c)] ;
n_var = numel(x) ;
x = reshape(x,n_var,1) ;
%% find empty and full elements
emptyelts = zeros(nely,nelx) ; fullelts = zeros(nely,nelx) ;
for i = 1:nely
    for j = 1:nelx
        %         % black round
        %         if sqrt((i-R)^2+(j-d)^2) <= R
        %             fullelts(i,j) = 1 ;
        %         end
        % white external round
        if sqrt((i-R)^2+(j-d)^2) > R
            emptyelts(i,j) = 1 ;
        end
    end
end
for i = 1:nely
    for j = 1:nelx
        % white internal round 1
        if sqrt((i-((3*nely)/5))^2+(j-(nelx/4))^2) <= Rint1
            emptyelts(i,j) = 1 ;
        end
    end
end
for i = 1:nely
    for j = 1:nelx
        % white internal round 2
        if sqrt((i-((3*nely)/5))^2+(j-(nelx/2))^2) <= Rint2
            emptyelts(i,j) = 1 ;
        end
    end
end
for i = 1:nely
    for j = 1:nelx
        % white internal round 3
        if sqrt((i-((3*nely)/5))^2+(j-((3*nelx)/4))^2) <= Rint3
            emptyelts(i,j) = 1 ;
        end
    end
end
for i = ((nely/10)+1):(9*nely)/10
    for j = f:g
        % white square
        emptyelts(i,j) = 1 ;
    end
end
emptyelts = find(emptyelts(:)) ;
fullelts = find(fullelts(:)) ;
% emptyelts = [] ; fullelts = [] ;
%% initialize plots
f1 = figure('units','normalized','position',[0.1,0.25,0.3,0.5]);
f2 = figure('units','normalized','position',[0.6,0.25,0.3,0.5]);
f3 = figure() ;
%% prepare finite element analysis
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
%% define loads and boundary conditions
% above
dinx1 = 2* (                    nely-a +1) -1;
diny1 = 2* (                    nely-a +1)   ;
dinx2 = 2* ((nely+1) * (d)             +1) -1;
diny2 = 2* ((nely+1) * (d)             +1)   ;
dinx3 = 2* ((nely+1) * (2*d)  + nely-a +1) -1;
diny3 = 2* ((nely+1) * (2*d)  + nely-a +1)   ;
dinx4 = 2* ((nely+1) * (nelx) + nely-b +1) -1;
diny4 = 2* ((nely+1) * (nelx) + nely-b +1)   ;
dout1 = 2* ((nely+1) * f   + (nely/10) +1)   ;
dout3 = 2* ((nely+1) * g   + (nely/10) +1)   ;
% below
dinx5 = 2* ((nely+1)            ) -1;
diny5 = 2* ((nely+1)            )   ;
dinx6 = 2* ((nely+1) * (d+1)    ) -1;
diny6 = 2* ((nely+1) * (d+1)    )   ;
dinx7 = 2* ((nely+1) * (2*d+1)  ) -1;
diny7 = 2* ((nely+1) * (2*d+1)  )   ;
dinx8 = 2* ((nely+1) * (nelx+1) ) -1;
diny8 = 2* ((nely+1) * (nelx+1) )   ;
dout2 = 2* ((nely+1) * f + (9*nely)/10);
dout4 = 2* ((nely+1) * g + (9*nely)/10);
% above
F(dinx1,1) =  1;
F(diny1,1) =  1;
F(dinx2,1) =  1;
F(diny2,1) =  1;
F(dinx3,1) =  1;
F(diny3,1) =  1;
F(dinx4,1) =  1;
F(diny4,1) =  1;
F(dout1,2) =  1;
F(dout3,2) =  1;
% below
F(dinx5,1) = -1;
F(diny5,1) = -1;
F(dinx6,1) = -1;
F(diny6,1) = -1;
F(dinx7,1) = -1;
F(diny7,1) = -1;
F(dinx8,1) = -1;
F(diny8,1) = -1;
F(dout2,2) = -1;
F(dout4,2) = -1;
% fixed/free degrees of freedom
fixeddofs1 = (2*(nely-a+1+1)-1):(2*(nely+1-1)-1);
fixeddofs2 = 2*(nely-a+1+1):2*(nely+1-1);
%fixeddofs1 = 2*(nely+1):2*(nely+1):2*(nely+1)*(nelx+1);
%fixeddofs2 = (2*(nely+1)-1):2*(nely+1):(2*(nely+1)*(nelx+1)-1);
fixeddofs = union(fixeddofs1,fixeddofs2);
alldofs     = (1:2*(nely+1)*(nelx+1));
freedofs    = setdiff(alldofs,fixeddofs);
U = zeros(2*(nely+1)*(nelx+1),2) ;
%% normalized variables
Xmin = repmat([0;0;-pi;0;0],n_c,1) ;
Xmax = repmat([nelx;nely;pi;nelx;nely],n_c,1) ;
x = (x-Xmin)./(Xmax-Xmin) ;
%% saving pictures
% parent_dir_name = 'mesh160x40_v0.45_16x4\';
% mkdir(parent_dir_name);
%% initialize optimization process
change = 1;
mm = 1;
n_var = length(x(:));
eeen    = ones(n_var,1);
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
kkttol  =tol;
kktnorm = kkttol+10;
loop = 0;
X_old=zeros(size(x));
X=Xmin+(Xmax-Xmin).*x;
%% start optimization process
while change > tol && loop <= maxoutit && kktnorm > tol
    loop = loop+1 ;
    x_old = x ;
    % density and it's gradient calculation
    den = zeros(nelx*nely,1); % density
    grad_den = zeros(nelx*nely,n_var); % gradient of density
    for in=1:n_c  % loop on the components
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
            den(ind) = den(ind) + f;
            grad_den(ind,n) = - cn*dfdx + sn*dfdy;
            grad_den(ind,n+1) = - sn*dfdx - cn*dfdy;
            grad_den(ind,n+2) = (-(xN(nc)-X(n))*sn + (yN(nc)-X(n+1))*cn)'.*dfdx - ...
                ((xN(nc)-X(n))*cn + (yN(nc)-X(n+1))*sn)'.*dfdy;
            grad_den(ind,n+3) = - ((xN(nc)-X(n))*cn + (yN(nc)-X(n+1))*sn)'.*dfdx/X(n+3);
            grad_den(ind,n+4) = ((xN(nc)-X(n))*sn - (yN(nc)-X(n+1))*cn)'.*dfdy/X(n+4);
        end
    end
    % ASYMPTOTIC DENSITY
    grad_den(den<1,:) = diag((1+12*den(den<1).^2-28*den(den<1).^3+15*den(den<1).^4))*grad_den(den<1,:);
    grad_den(den>=1,:) = 0;
    den = den+4*den.^3-7*den.^4+3*den.^5;
    den(fullelts) = 1; grad_den(fullelts,:) = 0;
    den(emptyelts) = 0; grad_den(emptyelts,:) = 0;
    den = sparse(min(den,1));
    grad_den = sparse(grad_den);
    % FEM analysis
    sK = reshape(KE(:)*(Emin+den'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
    % compliance and it's gradient calculation
    U1 = U(:,1) ; U2 = U(:,2) ;
    ce1 = sum((U1(edofMat)*KE).*U1(edofMat),2);
    ce2 = sum((U2(edofMat)*KE).*U2(edofMat),2);
    ce = ce1 + ce2 ;
    c = sum(sum((Emin+den.^penal*(E0-Emin)).*ce));
    dc = -penal*(E0-Emin)*grad_den'*(den.^(penal-1).*ce);
    cst=(sum(den(:))-mMax)/mMax*100;
    dcst=sum(grad_den)/mMax*100;
    dc=dc.*(Xmax-Xmin);dcst=dcst.*((Xmax-Xmin)');
    change=max(abs(X_old-X));
    X_old=X;
    f0val=c/10000;
    fval=cst;
    df0dx=dc(:)/10000;
    dfdx_=dcst(:)';
    % MMA code optimization
    [x(:),ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(mm,n_var,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx_,low,upp,a0,aa,C,d_);
    xold2 = xold1;
    xold1 = xval;
    xval  = x(:);
    fprintf(' It.:%5i Obj.:%11.4f Const.:%7.3f kktnorm.:%7.3f ch. %7.3f\n',loop,c, ...
        cst,kktnorm,change);
    % The residual vector of the KKT conditions is calculated:
    [~,kktnorm,~] = ...
        kktcheck(mm,n_var,x(:),ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx_,a0,aa,C,d_);
    % PLOT MEMBERS
    X=Xmin+(Xmax-Xmin).*x;
    figure(f1)
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
%     FileName=[parent_dir_name,'\Fig1_',int2str(loop),'.png'];
%     saveas(2,FileName);
    % plot densities
    figure(f2)
    colormap(flipud(gray)); imagesc(reshape(den,nely,nelx)); caxis([0 1]); axis equal; axis off; drawnow;
    title(['density field',int2str(loop)]) ;
    % convergence history
    if loop>0
        figure(f3)
        subplot(3,1,1)
        scatter(loop,c,4) ; xlabel('iteration') ; ylabel('compliance');
        title('compliance history') ;
        hold on
        subplot(3,1,2)
        scatter(loop,kktnorm,4) ; xlabel('iteration') ; ylabel('KKT norm');
        title('KKT norm history') ;
        hold on
        subplot(3,1,3)
        scatter(loop,cst,4) ; xlabel('iteration') ; ylabel('mass constraint (in %)');
        title('mass constraint history') ;
        hold on
    end
end