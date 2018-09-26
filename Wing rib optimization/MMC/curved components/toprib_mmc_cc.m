%% ************ Wing rib topology optimization with MNA **************** %%
% this code presents a topology optimization of a wing rib under MMC
% framework with curved components
clear ; close all ;
%% data
nelx = 160 ;
nely = 40 ;
DH = 1 ; DW = 4 ; EH = DH/nely ; EW = DW/nelx ;
volfrac = 0.45 ;
tol = 0.01 ;
E0 = 1 ; h = 1 ; Emin = 1e-6 ; penal = 3 ;
nu = 0.3 ;
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
l_ini = sqrt(volfrac*DW*DH/(0.5*nX*nY)) ; % initial length of component
[xElts,yElts] = meshgrid(linspace(1/(nX/2+1)*DW,nX/2/(nX/2+1)*DW,nX/2),...
    linspace(1/(nY/2+1)*DH,nY/2/(nY/2+1)*DH,nY/2)) ;
% removing components
for i = 1:length(xElts)
    % white external round
    if sqrt((yElts(i)-R*EH)^2+(xElts(i)-d*EH)^2) > R*EH
        xElts(i) = nan ; yElts(i) = nan ;
    end
end
for i = 1:length(xElts)
    % white internal round 1
    if sqrt((yElts(i)-((3*nely*EH)/5))^2+(xElts(i)-(nelx*EH/4))^2) <= Rint1*EH
        xElts(i) = nan ; yElts(i) = nan ;
    end
end
for i = 1:length(xElts)
    % white internal round 2
    if sqrt((yElts(i)-((3*nely*EH)/5))^2+(xElts(i)-(nelx*EH/2))^2) <= Rint2*EH
        xElts(i) = nan ; yElts(i) = nan ;
    end
end
for i = 1:length(xElts)
    % white internal round 3
    if sqrt((yElts(i)-((3*nely*EH)/5))^2+(xElts(i)-((3*nelx*EH)/4))^2) <= Rint3*EH
        xElts(i) = nan ; yElts(i) = nan ;
    end
end
for i = 1:length(xElts)
    % white square
    if f*EH<xElts(i) && xElts(i)<g*EH && ((nely*EH/10)+1)<yElts(i) && yElts(i)<(9*nely*EH)/10
        xElts(i) = nan ; yElts(i) = nan ;
    end
end
xElts = xElts(~isnan(xElts)) ;
yElts = yElts(~isnan(yElts)) ;
xElts = reshape(xElts,1,numel(xElts)) ;
yElts = reshape(yElts,1,numel(yElts)) ;
xElts = kron(xElts,[1 1]) ; yElts = kron(yElts,[1 1]) ;
n_c = length(xElts) ;
x = [xElts;yElts;repmat([pi/4 -pi/4],1,n_c/2);3*l_ini*ones(1,n_c);...
    1/3*l_ini*ones(1,n_c);0.2*ones(1,n_c)] ;
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
emptyelts = flipud(emptyelts) ;
emptyelts = find(emptyelts(:)) ;
fullelts = find(fullelts(:)) ;
% emptyelts = [] ; fullelts = [] ;
%% initialize plots
f1 = figure() ;
f2 = figure() ;
f3 = figure() ;
%% prepare finite element analysis
[KE] = BasicKe(E0,nu, EW, EH,h);
nodenrs=reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec=reshape(2*nodenrs(1:end-1,1:end-1)-1,nelx*nely,1);
edofMat=repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 4 5] 2 3],nelx*nely,1);
iK=kron(edofMat,ones(8,1))';
jK=kron(edofMat,ones(1,8))';
EleNodesID = edofMat(:,2:2:8)./2 ;
iEner = EleNodesID' ;
emptynodes = EleNodesID(emptyelts,:) ;
emptynodes = emptynodes(:) ;
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
Xmin = repmat([0;0;-pi;0;0;-1],n_c,1) ;
Xmax = repmat([DW;DH;pi;DW;DH;1],n_c,1) ;
x = (x-Xmin)./(Xmax-Xmin) ;
%% initialize optimization process
change = 1 ;
mm = 1 ;
n_var = length(x(:)) ;
eeen    = ones(n_var,1) ;
eeem    = ones(mm,1) ;
zerom   = zeros(mm,1) ;
xval    = x(:) ;
xold1   = xval ;
xold2   = xval ;
xmin    = 0*eeen ;
xmax    = eeen ;
low     = xmin ;
upp     = xmax ;
C       = 1000*eeem ;
d_       = 0*eeem ;
a0      = 1 ;
aa       = zerom ;
outeriter = 0;
maxoutit  = 300;
kkttol  =tol;
kktnorm = kkttol+10;
loop = 0;
X = Xmin+(Xmax-Xmin).*x ;
Phi = cell(n_c,1) ; % TDF of components
[mesh.x,mesh.y] = meshgrid(EW*(0:nelx),EH*(0:nely)) ; % meshing
p = 6 ; % parameter of TDF
alpha = 1e-2 ; % void parameter in the Heaviside function
epsilon = 4*min(EW,EH) ; % regularization parameter in Heaviside function
% epsilon = 4 ;
%% start optimization process
while change > tol && loop <= maxoutit && kktnorm > tol
    if loop > 15
        alpha = 1e-3 ;
    end
    loop = loop+1 ;
    X_old = X ;
    x_old = x ;
    % computing TDF for each component
    for i = 1:n_c
        Phi{i} = TDF_cc(X(6*(i-1)+1:6*i),mesh.x(:),mesh.y(:),p) ;
    end
    % computing global TDF
    P = 20 ; % parameter for KS function
    Phi_KS = zeros(size(Phi{1})) ;
    tempPhi_max = Phi{1} ;
    for i = 2:n_c
        tempPhi_max = max(tempPhi_max,Phi{i}) ;
    end
    for i = 1:n_c
        Phi_KS = Phi_KS+1/n_c*exp(P*(Phi{i}-tempPhi_max)) ;
    end
    Phi_KS = tempPhi_max+1/P*log(Phi_KS) ;
    Phi_KS(emptynodes) = -1 ;
    Phi_max = reshape(Phi_KS,nely+1,nelx+1) ;
    % plot components
    figure(f1) ;
    contourf(mesh.x,mesh.y,Phi_max,[0,0]) ;
    axis equal; axis([0 DW 0 DH]);  drawnow; title('component plot') ;
    % save pictures
    %     FileName=[parent_dir_name,'\Fig1_',int2str(Loop),'.png'];
    %     saveas(h,FileName);
    % Calculating the finite difference quotient of H
    H = Heaviside(Phi_max,alpha,nelx,nely,epsilon) ;
    H(emptynodes) = alpha ;
    diffH = cell(n_c,1) ;
    delta = max(2*min(EW,EH),0.00001) ;
%     delta = 0.1 ;
    for j = 1:n_c
        for ii = 1:6
            x1 = x ;
            x1(ii+(j-1)*6) = x(ii+(j-1)*6) + delta ;
            X1 = Xmin+x1.*(Xmax-Xmin) ;
            Phi1 = TDF_cc(X1((j-1)*6+1:6*j),mesh.x(:),mesh.y(:),p) ;
            Phi_max1 = Phi1 ;
            for s = 1:j-1
                Phi_max1 = max(Phi_max1,Phi{s}) ;
            end
            for s = j+1:n_c
                Phi_max1 = max(Phi_max1,Phi{s}) ;
            end
            Phi_KS1 = zeros(size(Phi{1})) ;
            for s = 1:j-1
                Phi_KS1 = Phi_KS1+1/n_c*exp(P*(Phi{s}-Phi_max1)) ;
            end
            for s = j+1:n_c
                Phi_KS1 = Phi_KS1+1/n_c*exp(P*(Phi{s}-Phi_max1)) ;
            end
            Phi_KS1 = Phi_KS1+1/n_c*exp(P*(Phi1-Phi_max1)) ;
            Phi_KS1 = Phi_max1+1/P*log(Phi_KS1) ;
            x2 = x ;
            x2(ii+(j-1)*6) = x(ii+(j-1)*6) - delta ;
            X2 = Xmin+x2.*(Xmax-Xmin) ;
            Phi2 = TDF_cc(X2((j-1)*6+1:6*j),mesh.x(:),mesh.y(:),p) ;
            Phi_max2 = Phi2 ;
            for s = 1:j-1
                Phi_max2 = max(Phi_max2,Phi{s}) ;
            end
            for s = j+1:n_c
                Phi_max2 = max(Phi_max2,Phi{s}) ;
            end
            Phi_KS2 = zeros(size(Phi{1})) ;
            for s = 1:j-1
                Phi_KS2 = Phi_KS2+1/n_c*exp(P*(Phi{s}-Phi_max2)) ;
            end
            for s = j+1:n_c
                Phi_KS2 = Phi_KS2+1/n_c*exp(P*(Phi{s}-Phi_max2)) ;
            end
            Phi_KS2 = Phi_KS2+1/n_c*exp(P*(Phi2-Phi_max2)) ;
            Phi_KS2 = Phi_max2+1/P*log(Phi_KS2) ;
            H1 = Heaviside(Phi_KS1,alpha,nelx,nely,epsilon) ;
            H2 = Heaviside(Phi_KS2, alpha,nelx,nely,epsilon) ;
            diffH{j}(:,ii) = (H1-H2)/(2*delta) ;
%             if ii == 1 || ii == 2 || ii == 3
%                 diffH{j}(:,ii) = 2*diffH{j}(:,ii) ;
%             else
%                 diffH{j}(:,ii) = 0.5*diffH{j}(:,ii) ;
%             end
            diffH{j}(emptynodes,:) = 0 ;
        end
    end
    % finite elements analysis
    denk = sum( H(EleNodesID).^2, 2 ) / 4;
    den=sum( H(EleNodesID), 2 ) / 4;
    grad_den = zeros(nelx*nely,n_var) ;
    for j = 1:n_c
        for i = 1:6
            dH_dxi = diffH{j}(:,i) ;
            grad_den(:,6*(j-1)+i) = 1/4*sum(dH_dxi(EleNodesID),2) ;
        end
    end
    A1=sum(den)*EW*EH ;
    sK = reshape(KE(:)*(Emin+den'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK(:),jK(:),sK(:)); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    %Energy of elememt and compliance
    U1 = U(:,1) ; U2 = U(:,2) ;
    % compliance and it's gradient calculation
    ce1 = sum((U1(edofMat)*KE).*U1(edofMat),2);
    ce2 = sum((U2(edofMat)*KE).*U2(edofMat),2);
    ce = ce1 + ce2 ;
    c = sum(sum((Emin+den.^penal*(E0-Emin)).*ce));
    dc = -penal*(E0-Emin)*grad_den'*(den.^(penal-1).*ce);
    % sensitivities
    dfdx = zeros(n_var,1) ;
    for k = 1:n_c
        dfdx(6*(k-1)+1:6*k) = sum(diffH{k})/4 ;
    end
    f0val = c/10000 ;
    df0dx = dc/10000 ;
    fval = (A1/(DW*DH)-volfrac)/volfrac ;
    dfdx = dfdx*EW*EH/(DW*DH)/volfrac ;
    % MMA code optimization
    [x(:),ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(mm,n_var,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,aa,C,d_);
    xold2 = xold1 ;
    xold1 = xval ;
    xval  = x(:) ;
    X = Xmin + (Xmax-Xmin).*xval ;
    change = max(abs(X-X_old)) ;
    % The residual vector of the KKT conditions is calculated:
    [~,kktnorm,~] = ...
        kktcheck(mm,n_var,x(:),ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,aa,C,d_);
    % display iteration's results
    fprintf(' It.:%5i Obj.:%11.4f Const.:%7.3f kktnorm.:%7.3f ch. %7.3f\n',loop,c, ...
        fval,kktnorm,change) ;
    % plot densities
    figure(f2)
    colormap(flipud(gray)); imagesc(flipud(reshape(den,nely,nelx))); caxis([0 1]); axis equal; axis off; drawnow;
    title('density field') ;
    % convergence history
    if loop>20
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
        scatter(loop,fval,4) ; xlabel('iteration') ; ylabel('mass constraint (in %)');
        title('mass constraint history') ;
        hold on
    end
end