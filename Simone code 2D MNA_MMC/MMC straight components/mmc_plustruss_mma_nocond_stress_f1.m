%% ***************** truss model code with MMC method ****************** %%
clear ; close all ; tic
%% DATA
DH = 1 ;
DW = 2 ;
nelx = 80 ;
nely = 40 ;
penal = 3 ;
x_int = 0.25 ;
y_int = 0.25 ;
ini_val = [0.38 0.04 0.06 0.04 0.7] ;
C0 = 70 ;
T0 = 15 ;
M = [ nely + 1 , nelx + 1 ];
EW = DW / nelx;  % length of element
EH = DH / nely;  % width of element
[ x ,y ] = meshgrid( EW * [ 0 :  nelx] , EH * [ 0 : nely]); %#ok<*NBRAK>
LSgrid.x = x(:);
LSgrid.y = y(:);     % coordinate of nodes
%% MATERIAL PROPERTIES
E = 1;
nu = 0.3;
P=4;
Sl=1;
%% Component geometry initialization
x0=x_int/2:x_int:DW;    % x-coordinates of the centers of components
y0=y_int/2:y_int:DH;    % y-coordinates of the centers of components
xn=length(x0);               % number of component groups in x direction
yn=length(y0);               % number of component groups in y direction
x0=kron(x0,ones(1,2*yn));
y0=repmat(kron(y0,ones(1,2)),1,xn);
N=length(x0);                % total number of components in the design domain
L=repmat(ini_val(1),1,N);                    % vector of the half length of each component
t1=repmat(ini_val(2),1,N);                   % vector of the half width of component at point A
t2=repmat(ini_val(3),1,N);                   % vector of the half width of component at point B
t3=repmat(ini_val(4),1,N);                   % vector of the half width of component at point C
st=repmat([ini_val(5) -ini_val(5)],1,N/2);   % vector of the sine value of the inclined angle of each component
variable=[x0;y0;L;t1;t2;t3;st];
%% Parameters of MMA
xy=variable(:);
xval=xy;
xold1 = xy;
xold2 = xy;
%% PREPARE FINITE ELEMENT ANALYSIS
[KE] = BasicKe(E,nu, EW, EH,1);   %  stiffness matrix k^s is formed
D=E/(1-nu^2)*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
B=1/2*[-1 0 1 0 1 0 -1 0;0 -1 0 -1 0 1 0 1;-1 -1 -1 1 1 1 1 -1];
DB=D*B;
Cvm=[1 -0.5 0;-0.5 1 0;0 0 3];
Sel=DB'*Cvm*DB;
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
EleNodesID=edofMat(:,2:2:8)./2;
iEner=EleNodesID';
fixeddofs = [2*(fix(nelx/2)*(nely+1)+1:(nely+1):((nelx)*(nely+1)+1))-1;2*(fix(nelx/2)*(nely+1)+1:(nely+1):((nelx)*(nely+1)+1))];
fixeddofs=fixeddofs(:);

%% Generate the truss stiffness matrix
A=1;
I=100;
[Kcc,Kce,Kee,Fe,Fc,Recovery_matrixc,Recovery_matrixe]=truss_stiffness_no_condensation(nelx,nely,E,A,I);
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
%% saving pictures
parent_dir_name = 'mesh80x40\';
mkdir(parent_dir_name);
%% INITIALIZE ITERATION
Xmin    = repmat([0;0;0.01;0.01;0.01;0.03;-1],length(xy)/7,1);
Xmax    = repmat([DW;DH;2;0.2;0.2;0.2;1],length(xy)/7,1);
xy = (xy-Xmin)./(Xmax-Xmin) ;
Var_num=7;  % number of design variables for each component
m = 3;
n = length(xy);
eeen    = ones(n,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xmin    = zeron;
xmax    = eeen;
low     = xmin;
upp     = xmax;
C       = 1000*eeem;
D       = 0*eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 2300;
kkttol  =0.01;
p = 6 ;
alpha=1e-2;            % parameter alpha in the Heaviside function
epsilon=4*min(EW,EH);  % regularization parameter in the Heaviside function
kktnorm = kkttol+10;
outit = 0;
[ud,Sd,~]=svds(Sel,8);
X_old=zeros(size(xy));
X=Xmin+(Xmax-Xmin).*xy;
Phi=cell(N,1);
%% START ITERATION
while kktnorm > kkttol && outit < maxoutit
    outit   = outit+1;
    outeriter = outeriter+1;
    %Forming Phi^s
    for i=1:N
        Phi{i}=tPhi(X(Var_num*i-Var_num+1:Var_num*i),LSgrid.x,LSgrid.y,p);
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
    Phi_max=reshape(Phi_KS,nely+1,nelx+1);
    %Plot components
    figure (2)
    contourf(reshape(x , M), reshape(y , M),flipud(Phi_max),[0,0]);
    axis equal;axis([0 DW 0 DH]);
    title(['contour plot at iteration ',num2str(outit)]);pause(1e-6); 
    % save pictures
    FileName=[parent_dir_name,'\Fig1_',int2str(outit),'.png'];
    saveas(1,FileName);
    % Calculating the finite difference quotient of H
    H=Heaviside(Phi_max,alpha,nelx,nely,epsilon);
    diffH=cell(N,1);
    delta=max(2*min(EW,EH),0.00001);
    for j=1:N
        for ii=1:Var_num
            X001=X;
            X001(ii+(j-1)*Var_num)=X(ii+(j-1)*Var_num)+delta;
            tmpPhiD1=tPhi(X001(Var_num*j-Var_num+1:Var_num*j),LSgrid.x,LSgrid.y,p);
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
            X002=X;
            X002(ii+(j-1)*Var_num)=X(ii+(j-1)*Var_num)-delta;
            tmpPhiD2=tPhi(X002(Var_num*j-Var_num+1:Var_num*j),LSgrid.x,LSgrid.y,p);
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
        end
    end
    %% FE-ANALYSIS
    denk = sum( H(EleNodesID).^2, 2 ) / 4;
    den=sum( H(EleNodesID), 2 ) / 4;
    % density gradient
    den_grad = zeros(nelx*nely,n) ;
    for j=1:N
        for ii=1:Var_num
            dH = diffH{j}(:,ii) ;
            den_grad(:,ii+(j-1)*Var_num) = sum(dH(EleNodesID),2)/4 ;
        end
    end
    A1=sum(den)*EW*EH;
    sK = KE(:)*denk(:)';
    K0 = sparse(iK,jK,sK,size(Kt,1),size(Kt,2)); K=K0+Kt;K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    %Energy of elememt
    energy = sum((U(edofMat)*KE).*U(edofMat),2);
    sEner=ones(4,1)*energy'/4;
    energy_nod=sparse(iEner(:),1,sEner(:));
    %total compliance
    compliance=F'*U;
    % Sensitivities
    dcompliance=zeros(Var_num*N,1);
    df0dx=zeros(Var_num*N,1);
    for k=1:N
        dcompliance(Var_num*k-Var_num+1:Var_num*k,1)=2*energy_nod'.*H*diffH{k};
        df0dx(Var_num*k-Var_num+1:Var_num*k,1)=sum(diffH{k})/4;
    end
    %relative squared displacement
    mse = reshape(sqrt(sum(((U(edofMat)*ud(:,1:3))*Sd(1:3,1:3)).*(U(edofMat)*ud(:,1:3)),2)),nely,nelx); %microscopic Von Mises Stress
    rcv=reshape(den,nely,nelx).*(mse/Sl-1);
    G_ks=max(rcv(:))+1/P*log(mean(exp(P*(rcv(:)-max(rcv(:))))));
    rcv=den.*(mse(:)/Sl-1);
    sS=reshape(Sel(:)*(den'./mse(:)'/Sl.*exp(P*(rcv(:)'-max(rcv(:)))))/sum(exp(P*(rcv(:)-max(rcv(:))))),64*nelx*nely,1);
    S0=sparse(iK,jK,sS,length(F),length(F)); S0=(S0+S0')/2;
    dG_ksdu=S0*U;
    rd=(Rt*U);
    c=((rd')*rd)/length(rd);
    Lambda(freedofs,:)=K(freedofs,freedofs)\[(2*Rt(:,freedofs)'*rd),dG_ksdu(freedofs)];
    Lambda1=Lambda(:,1);
    Lambda2=Lambda(:,2);
    dG_ksdx=den_grad'*((mse(:)/Sl-1).*exp(P*(rcv(:)-max(rcv(:))))/sum(exp(P*(rcv(:)-max(rcv(:))))));
    dG_du = sum((U(edofMat)*KE).*Lambda2(edofMat),2) ;
    sdG_du=ones(4,1)*dG_du'/4;
    dG_du_nod=sparse(iEner(:),1,sdG_du(:));
    ce = sum((Lambda1(edofMat)*KE).*U(edofMat),2) ;
    sce=ones(4,1)*dG_du'/4;
    ce_nod=sparse(iEner(:),1,sce(:));
    dG_du = zeros(n,1) ;
    for k=1:N
        dG_du(Var_num*k-Var_num+1:Var_num*k,1)=2*dG_du_nod'.*H*diffH{k};
        dc(Var_num*k-Var_num+1:Var_num*k,1)=-2/length(rd)*ce_nod'.*H*diffH{k};
    end
    dGKSl=dG_du(:)+dG_ksdx ;
    %% Gather info for MMA
    f0val=mean(den)*100;
    fval=[(c-T0)/T0*100;G_ks*100;(compliance-C0)/C0*100];%;1*(-0.38-FAN_Ax_disp)/0.38
    fval=fval/1000;
    df0dx=mean(den_grad)'*100;
    dfdx=[dc(:)'/T0*100;dGKSl(:)'*100;dcompliance(:)'/C0*100];%;-1*dfandisp(:)'/0.38
    dfdx=-dfdx/1000;
    %% MMA code optimization
    [xy(:),ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,C,D);
    xold2 = xold1;
    xold1 = xval;
    xval  = xy(:);
    change=norm(xval-xold1);
    %% %% The residual vector of the KKT conditions is calculated:
    [~,kktnorm,~] = ...
        kktcheck(m,n,xy(:),ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,C,D);
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f kktnorm.:%7.3f ch.:%7.3f\n',outit,c, ...
        full(f0val),kktnorm,change);
    %% PLOT DENSITIES
    figure(1)
    h = figure(1); set(h,'Color',[1 1 1]);
    [Yy,Xx]=find(nodenrs);
    Yy=nely+1-Yy;
    Xx=Xx-1;
    hold on;  patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-den)*[1 1 1],'FaceColor','flat','EdgeColor','none'); axis equal; axis off; drawnow;hold off
    title(['Design zone at iteration ',num2str(outit)])
    colormap(gray);
    % non normalized variables
    X=Xmin+(Xmax-Xmin).*xy;
    % convergence history
    if outit>5
        figure(3)
        subplot(3,2,1)
        scatter(outit,f0val,4) ; xlabel('iteration') ; ylabel('volume');
        title('volume history') ;
        hold on
        subplot(3,2,2)
        scatter(outit,kktnorm,4) ; xlabel('iteration') ; ylabel('KKT norm');
        title('KKT norm history') ;
        hold on
        subplot(3,2,3)
        scatter(outit,1000*full(fval(1)),4) ; xlabel('iteration') ; ylabel('TSFC constraint (in %)');
        title('TSFC constraint history') ;
        hold on
        subplot(3,2,4)
        scatter(outit,10*full(fval(2)),4) ; xlabel('iteration') ; ylabel('aggregated stress constraint');
        title('stress constraint history') ;
        hold on
        subplot(3,2,5)
        scatter(outit,1000*full(fval(3)),4) ; xlabel('iteration') ; ylabel('compliance constraint (in %)');
        title('compliance constraint history') ;
        hold on
    end
end
toc