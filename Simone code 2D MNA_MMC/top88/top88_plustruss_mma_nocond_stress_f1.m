%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function [x,compliance,c,f0val,outit,G_ks]=top88_plustruss_mma_nocond_stress_f1(nelx,nely,penal,rmin,ft,C0,T0)
close all
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

%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

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
x = ones(nely,nelx);
xPhys = x;
m = 3;
n = length(xPhys(:));
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
C       = 1000*eeem;
d       = 0*eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 10000;
kkttol  =0.01;
%
%%%% The iterations start:
kktnorm = kkttol+10;
% kktnorm = kkttol;
outit = 0;
[ud,Sd,~]=svds(Sel,8);

%% START ITERATION
while kktnorm > kkttol && outit < maxoutit
    outit   = outit+1;
    outeriter = outeriter+1;
    if ft == 1
        xPhys(:)=xval;
    elseif ft == 2
        xPhys(:) = (H*xval(:))./Hs;
    end
    %% FE-ANALYSIS
    
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K0 = sparse(iK,jK,sK,size(Kt,1),size(Kt,2)); K=K0+Kt;K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
%      cndK=log(condest(K(freedofs,freedofs)))/log(10);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    %compliance of DZ only
    %   Tau=K0*U;
    %   Lambda(freedofs)=K(freedofs,freedofs)\Tau(freedofs);
    %   ce = 2*reshape(sum((U(edofMat)*KE).*Lambda(edofMat),2),nely,nelx)-reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    %   c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
    %   dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx);
    %total compliance
    compliance=F'*U;
    ce_comp=reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    dcompliance = -penal*(E0-Emin)*xPhys.^(penal-1).*ce_comp;
    %relative squared displacement
%     rd1=(Gamma*Rt*U);

  mse = reshape(sqrt(sum(((U(edofMat)*ud(:,1:3))*Sd(1:3,1:3)).*(U(edofMat)*ud(:,1:3)),2)),nely,nelx); %microscopic Von Mises Stress
  rcv=xPhys.*(mse/Sl-1);
  G_ks=max(rcv(:))+1/P*log(mean(exp(P*(rcv(:)-max(rcv(:))))));
  rcv=xPhys(:).*(mse(:)/Sl-1);
  sS=reshape(Sel(:)*(xPhys(:)'./mse(:)'/Sl.*exp(P*(rcv(:)'-max(rcv(:)))))/sum(exp(P*(rcv(:)-max(rcv(:))))),64*nelx*nely,1);
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
  dG_ksdx=(mse(:)/Sl-1).*exp(P*(rcv(:)-max(rcv(:))))/sum(exp(P*(rcv(:)-max(rcv(:)))));
  dG_du=reshape((sum((U(edofMat)*KE).*Lambda2(edofMat),2)),nely,nelx);
  dG_du =- penal*(E0-Emin)*xPhys(:).^(penal-1).*dG_du(:);
  dGKSl=dG_du(:)+dG_ksdx;
    ce=reshape(sum((Lambda1(edofMat)*KE).*U(edofMat),2),nely,nelx);
    dc = -1/length(rd)*penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
        dGKSl(:) = H*(dGKSl(:)./Hs);
    end
    %% Gather info for MMA
    f0val=mean(xPhys(:))*100;
    fval=[(c-T0)/T0*100;G_ks*100;(compliance-C0)/C0*100];%;1*(-0.38-FAN_Ax_disp)/0.38
    df0dx=dv(:)/length(dv(:))*100;
    dfdx=[dc(:)'/T0*100;dGKSl(:)'*100;dcompliance(:)'/C0*100];%;-1*dfandisp(:)'/0.38
    %% MMA code optimization
    [x(:),ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,C,d);
    xold2 = xold1;
    xold1 = xval;
    xval  = x(:);
    change=norm(xval-xold1);
    %     xval(x1>=2500&x4<=4500&z1==min(z1))=0;
    
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f kktnorm.:%7.3f ch.:%7.3f\n',outit,c, ...
        mean(xPhys(:)),kktnorm,change);
%     figure(2)
%     hold on
%     plot(outit,c,'bo','MarkerFaceColor','b')
%     plot(outit,mean(xPhys(:))*100,'ro','MarkerFaceColor','r')
%     plot(outit,G_ks*100,'mo','MarkerFaceColor','m')
%     
%     plot(outit,compliance,'go','MarkerFaceColor','g')
%     % plot(outeriter,(1+GKSl)*VMl,'ko','MarkerFaceColor','k')
%     title(['Convergence V = ',num2str(mean(xPhys(:))*100),', T =',num2str(c),', iter = ', num2str(outit)])
%     grid on
%     legend('T','V %',['G_{KS}^l = ',num2str(G_ks*100),'%'],['c =',num2str(compliance)])
%     xlabel('iter')
    %% %% The residual vector of the KKT conditions is calculated:
    [~,kktnorm,~] = ...
        kktcheck(m,n,x(:),ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,C,d);
    
    %% PLOT DENSITIES
      figure(1)
%       colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    h = figure(1); set(h,'Color',[1 1 1]);
    [Yy,Xx]=find(nodenrs);
    Yy=nely+1-Yy;
    Xx=Xx-1;
    hold on;  patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xPhys(:))*[1 1 1],'FaceColor','flat','EdgeColor','none'); axis equal; axis off; drawnow;hold off
    title(['Design zone at iteration ',num2str(outit)])
    colormap(gray);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


