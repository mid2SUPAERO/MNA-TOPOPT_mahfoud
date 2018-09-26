%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function x=top88_mma_multi_resol(nelx,nely,volfrac,penal,rmin,ft)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
nl=1;
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
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
ndof=2*(nelx+1)*(nely+1);
Pu=speye(ndof);
for s=1:nl
[Pu] = Pu*Prepcoarse(nely/2^(s-1),nelx/2^(s-1));
end
N = ones(ndof,1); N(fixeddofs) = 0; Null = spdiags(N,0,ndof,ndof);
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
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
xPhys = x;
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
C       = 100000*eeem;
d       = 0*eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 1500;
kkttol  =0.001;
%
%%%% The iterations start:
kktnorm = kkttol+10;
% kktnorm = kkttol;
outit = 0;
change=1;
Fc=Pu'*F;
%% START ITERATION
while kktnorm > kkttol && outit < maxoutit && change>0.01
    outit   = outit+1;
    outeriter = outeriter+1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  K = Null'*K*Null - (Null-speye(ndof,ndof));
  tic
  Kc=Pu'*K*Pu;
  eval_k=toc
  tic
  U = Pu*(Kc\Fc);
  eval_inv=toc
%   cndK=log(condest(K(freedofs,freedofs)))/log(10);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  dv = ones(nely,nelx);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  %% Gather info for MMA
f0val=c;
fval=mean(xPhys(:))-volfrac;%;1*(-0.38-FAN_Ax_disp)/0.38
df0dx=dc(:);
dfdx=dv(:)'/length(dv(:));%;-1*dfandisp(:)'/0.38
innerit=0;
outvector1 = [outeriter innerit f0val fval'];
outvector2 = xval;

  %% MMA code optimization
    [x(:),ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,C,d);
    xold2 = xold1;
    xold1 = xval;
    xval  = x(:);
    change=norm(xval-xold1);
    %     xval(x1>=2500&x4<=4500&z1==min(z1))=0;
    if ft == 1
        xPhys(:)=xval;
    elseif ft == 2
        xPhys(:) = (H*xval(:))./Hs;
    end
    %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f kktnorm.:%7.3f ch.:%7.3f\n',outit,full(c), ...
    mean(xPhys(:)),kktnorm,change);
figure(2)
hold on
plot(outit,c,'bo','MarkerFaceColor','b')
plot(outit,mean(xPhys(:))*100,'ro','MarkerFaceColor','r')
% plot(outeriter,(1+GKSl)*VMl,'ko','MarkerFaceColor','k')
title(['Convergence volfrac = ',num2str(mean(xPhys(:))*100),', Compliance =',num2str(gather(c)),', iter = ', num2str(outit)])
grid on
legend('Compliance','Volume Fraction %')
xlabel('iter')
 %% %% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(m,n,x(:),ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,C,d);
    outvector1 = [outeriter innerit f0val fval(:)'];
    outvector2 = xval;
  
  %% PLOT DENSITIES
  figure(1)
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
%
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


