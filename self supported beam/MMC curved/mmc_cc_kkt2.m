%% MMC code with curved components fixed width %%
clear all ; close all ; tic %#ok<*CLALL>
%% data
DH = 1 ;
DW = 2 ;
nelx = 160 ;
nely = 80 ;
x_int = 0.5 ;
y_int = 0.5 ;
ini_val = [0.4 0.04 pi/4 0] ;
volfrac = 0.5 ;
%% saving pictures
% parent_dir_name = 'mesh80x40_old_var_v0.5\';
% mkdir(parent_dir_name);
%% FEM data initialization
M = [ nely + 1 , nelx + 1 ];
EW = DW / nelx;  % length of element
EH = DH / nely;  % width of element
[ x ,y ] = meshgrid( EW * [ 0 :  nelx] , EH * [ 0 : nely]); %#ok<*NBRAK>
LSgrid.x = x(:);
LSgrid.y = y(:);     % coordinate of nodes
%% Material properties
h=1;  %thickness
E=1;
nu=0.3;
%% Component geometry initialization
x0=x_int/2:x_int:DW;    % x-coordinates of the centers of components
y0=y_int/2:y_int:DH;    % y-coordinates of the centers of components
xn=length(x0);               % number of component groups in x direction
yn=length(y0);               % number of component groups in y direction
x0=kron(x0,ones(1,2*yn));
y0=repmat(kron(y0,ones(1,2)),1,xn);
% x0=kron(x0,ones(1,yn));
% y0=repmat(kron(y0,ones(1,1)),1,xn);
N=length(x0);                % total number of components in the design domain
L=repmat(ini_val(1),1,N);                    % vector of the half length of each component
t=repmat(ini_val(2),1,N);                   % vector of the half width of component
theta=repmat([ini_val(3) -ini_val(3)],1,N/2); % vector of the sine value of the inclined angle of each component
% R=repmat([ini_val(4) -ini_val(4)],1,N/2);
% R=repmat(ini_val(4),1,N);
% cr=repmat(ini_val(4),1,N);
cr=repmat([ini_val(4) -ini_val(4)],1,N/2);
variable=[x0;y0;L;t;theta;cr];
N = length(x0) ;
%% Parameters of MMA
xy00=variable(:);
xval=xy00;
xold1 = xy00;
xold2 = xy00;
%% Limits of variable:[x0 ; y0  ; L ; t ; theta ; cr];
xmin=[0 ; 0  ; 0.03 ; 0.01 ; -2*pi ; -1];
xmin=repmat(xmin,N,1);
xmax=[DW ; DH ; 1 ; 0.2 ; 2*pi ; 1];
xmax=repmat(xmax,N,1);
low   = xmin;
upp   = xmax;
m = 1;  %number of constraint
Var_num=6;  % number of design variables for each component
nn=Var_num*N;
c=1000*ones(m,1);
d=zeros(m,1);
a0=1;
a=zeros(m,1);
%% Define loads and supports(self supported beam)
fixeddofs = [2 2*(nely+1)*nelx+2] ; 
alldofs     = 1:2*(nely+1)*(nelx+1) ; 
freedofs    = setdiff(alldofs,fixeddofs) ; 
loaddof = (nely+1)*nelx+nely ;
%% Preparation FE analysis
nodenrs=reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec=reshape(2*nodenrs(1:end-1,1:end-1)-1,nelx*nely,1);
edofMat=repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 4 5] 2 3],nelx*nely,1);
iK=kron(edofMat,ones(8,1))';
jK=kron(edofMat,ones(1,8))';
EleNodesID=edofMat(:,2:2:8)./2;
iEner=EleNodesID';
[KE] = BasicKe(E,nu, EW, EH,h);   %  stiffness matrix k^s is formed
% emptynodes = EleNodesID(emptyelts,:) ;
% emptynodes = emptynodes(:) ;
%% Initialize iteration
p=6;
alpha=1e-2;            % parameter alpha in the Heaviside function
epsilon=4*min(EW,EH);  % regularization parameter epsilon in the Heaviside function
Phi=cell(N,1);
Loop=1;
maxiter=600;      % the maximum number of iterations
kkttol=0.001;
kktnorm = kkttol+10;
list_comp=[];
list_vol=[];
list_kktnorm=[];
list_iter=[];
xy_opt=xy00; % initialisation du design optimal
while Loop<maxiter  && kktnorm>kkttol
    if Loop>20
        alpha = 1e-3 ;
    end
    %Forming Phi^s
    for i=1:N
        Phi{i}=tPhi_c(xy00(Var_num*i-Var_num+1:Var_num*i),LSgrid.x,LSgrid.y,p);
    end
    %Union of components KS approximation
    P=20; % parameter for KS function
    Phi_KS=zeros(size(Phi{1}));
    tempPhi_max=Phi{1};
    for i=2:N
        tempPhi_max=max(tempPhi_max,Phi{i});
    end
    for i=1:N
        Phi_KS=Phi_KS+1/N*exp(P*(Phi{i}-tempPhi_max));
    end
    Phi_KS=tempPhi_max+1/P*log(Phi_KS);
%     Phi_KS(emptynodes)=-1 ;
    Phi_max=reshape(Phi_KS,nely+1,nelx+1);
    %Plot components
    contourf(reshape(x , M), reshape(y , M),Phi_max,[0,0]);
    axis equal;axis([0 DW 0 DH]);pause(1e-6);
    % save pictures
    %     FileName=[parent_dir_name,'\Fig1_',int2str(Loop),'.png'];
    %     saveas(h,FileName);
    % Calculating the finite difference quotient of H
    H=Heaviside(Phi_max,alpha,nelx,nely,epsilon);
%     H(emptynodes)=alpha;
    diffH=cell(N,1);
    delta=max(2*min(EW,EH),0.00001);
%     delta=0.01;
    for j=1:N
        for ii=1:Var_num
            xy001=xy00;
            xy001(ii+(j-1)*Var_num)=xy00(ii+(j-1)*Var_num)+delta;
            tmpPhiD1=tPhi_c(xy001(Var_num*j-Var_num+1:Var_num*j),LSgrid.x,LSgrid.y,p);
            tempPhi_max1=tmpPhiD1;
            for ik=1:j-1
                tempPhi_max1=max(tempPhi_max1,Phi{ik});
            end
            for ik=j+1:N
                tempPhi_max1=max(tempPhi_max1,Phi{ik});
            end
            Phi_KS1=zeros(size(Phi{1}));
            for i=1:j-1
                Phi_KS1=Phi_KS1+1/N*exp(P*(Phi{i}-tempPhi_max1));
            end
            for i=j+1:N
                Phi_KS1=Phi_KS1+1/N*exp(P*(Phi{i}-tempPhi_max1));
            end
            Phi_KS1=Phi_KS1+1/N*exp(P*(tmpPhiD1-tempPhi_max1));
            Phi_KS1=tempPhi_max1+1/P*log(Phi_KS1);
            xy002=xy00;
            xy002(ii+(j-1)*Var_num)=xy00(ii+(j-1)*Var_num)-delta;
            tmpPhiD2=tPhi_c(xy002(Var_num*j-Var_num+1:Var_num*j),LSgrid.x,LSgrid.y,p);
            tempPhi_max2=tmpPhiD2;
            for ik=1:j-1
                tempPhi_max2=max(tempPhi_max2,Phi{ik});
            end
            for ik=j+1:N
                tempPhi_max2=max(tempPhi_max2,Phi{ik});
            end
            Phi_KS2=zeros(size(Phi{1}));
            for i=1:j-1
                Phi_KS2=Phi_KS2+1/N*exp(P*(Phi{i}-tempPhi_max2));
            end
            for i=j+1:N
                Phi_KS2=Phi_KS2+1/N*exp(P*(Phi{i}-tempPhi_max2));
            end
            Phi_KS2=Phi_KS2+1/N*exp(P*(tmpPhiD2-tempPhi_max2));
            Phi_KS2=tempPhi_max2+1/P*log(Phi_KS2);
            HD1=Heaviside(Phi_KS1,alpha,nelx,nely,epsilon);
            HD2=Heaviside(Phi_KS2,alpha,nelx,nely,epsilon);
            diffH{j}(:,ii)=(HD1-HD2)/(2*delta);
%             diffH{j}(emptynodes,:)=0;
        end
    end
    %FEA
    denk = sum( H(EleNodesID).^2, 2 ) / 4;
    den=sum( H(EleNodesID), 2 ) / 4;
    A1=sum(den)*EW*EH;
    F=sparse(loaddof,1,-A1/(DH*DW),2*(nely+1)*(nelx+1),1);
    U = zeros(2*(nely+1)*(nelx+1),1);
    sK = KE(:)*denk(:)';
    K = sparse(iK(:),jK(:),sK(:)); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    %Energy of elememt
    energy = sum((U(edofMat)*KE).*U(edofMat),2);
    sEner=ones(4,1)*energy'/4;
    energy_nod=sparse(iEner(:),1,sEner(:));
    Comp=F'*U;
    % Sensitivities
    dden = zeros(nelx*nely,nn) ;
    for j = 1:nn
        k = 1+floor((j-1)/Var_num) ; s = 1+rem((j-1),Var_num) ;
        dH_dx = diffH{k}(:,s) ;
        dden(:,j) = sum(dH_dx(EleNodesID),2)/4 ;
    end
    dg = -(EW*EH)/(DW*DH)*dden'*ones(nelx*nely,1) ;
    dF = zeros(2*(nelx+1)*(nely+1),nn) ;
    dF((nely+1)*nelx+nely,:) = dg' ; dF = sparse(dF) ;
    df0dx=zeros(Var_num*N,1);
    dfdx = sum(dden)/(nelx*nely*volfrac)*100 ; dfdx = dfdx' ;
%     dfdx=zeros(Var_num*N,1);
    for k=1:N
        df0dx(Var_num*k-Var_num+1:Var_num*k,1)=2*energy_nod'.*H*diffH{k};
%         dfdx(Var_num*k-Var_num+1:Var_num*k,1)=sum(diffH{k})/4;
    end
    %MMA optimization
    f0val = Comp/1000 ;
    %     df0dx=-df0dx/max(abs(df0dx));
    df0dx = (-df0dx + 2*dF'*U)/1000 ;
    fval = (A1/(DW*DH)-volfrac)/volfrac*100 ;
    %     dfdx=dfdx/max(abs(dfdx));
%     dfdx = -dfdx*EW*EH/(DW*DH)/volfrac*100 ;
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,ss,low,upp] = ...
        mmasub(m,nn,Loop,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d); %#ok<*ASGLU>
    xold2 = xold1;
    xold1 = xval;
    % The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(m,nn,xmma,ymma,zmma,lam,xsi,eta,mu,zet,ss, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
    %     outvector1 = [Loop Loop f0val fval(:)'];
    %     outvector2 = xval;
    xval = xmma;
    xy00=round(xval*1e4)/1e4;
    %     xy00=xval;
    % mise à jour des listes
    list_iter(Loop)=Loop; %#ok<*SAGROW>
    list_comp(Loop)=Comp;
    list_vol(Loop)=fval;
    list_kktnorm(Loop)=kktnorm;
    % mise à jour de la solution optimale
    if fval<=0 && Loop>20 && Comp<list_comp(Loop-1)
        xy_opt=xy00;
    end
    % affichage
    disp([' It.: ' sprintf('%4i\t',Loop) ' Obj.: ' sprintf('%6.3f\t',Comp) ' Vol.: ' ...
        sprintf('%6.4f\t',fval) 'kkt.:' sprintf('%6.4f\t',kktnorm)]);
    Loop = Loop + 1;
end
toc
%% historique de convergence
figure;
subplot(2,2,1)
semilogy(list_iter(30:end),list_comp(30:end))
xlabel('iteration') ; ylabel('compliance')
title('compliance history')
subplot(2,2,2)
plot(list_iter(30:end),list_vol(30:end))
xlabel('iteration') ; ylabel('volume')
title('volume history')
subplot(2,2,3)
semilogy(list_iter(30:end),list_kktnorm(30:end))
xlabel('iteration') ; ylabel('kktnorm')
title('kkt-norm history')
% solution optimale (min des itérations)
Phi_opt=cell(1,N);
%Forming Phi^s
for i=1:N
    Phi_opt{i}=tPhi_c(xy_opt(Var_num*i-Var_num+1:Var_num*i),LSgrid.x,LSgrid.y,p);
end
%Union of components KS approximation
Phi_KS_opt=zeros(size(Phi_opt{1}));
tempPhi_max_opt=Phi_opt{1};
for i=2:N
    tempPhi_max_opt=max(tempPhi_max_opt,Phi_opt{i});
end
for i=1:N
    Phi_KS_opt=Phi_KS_opt+1/N*exp(P*(Phi_opt{i}-tempPhi_max_opt));
end
Phi_KS_opt=tempPhi_max_opt+1/P*log(Phi_KS_opt);
% Phi_KS_opt(emptynodes)=-1 ;
Phi_max_opt=reshape(Phi_KS_opt,nely+1,nelx+1);
%Plot components
figure;
contourf(reshape(x , M), reshape(y , M),Phi_max_opt,[0,0]);
axis equal;axis([0 DW 0 DH]);