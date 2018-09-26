clear
close all
load('substructure_results')
load('st3D.dat.mat')
load('dz3Dstr.mat')
load('connectivity_interface')
load('gamma_structure')
load('U_tilder')
%% gather all needed input in FEM_structure
FEM_structure.DOFs_MAP=DOFs_MAP;
FEM_structure.K_interface=sparse(K_interface);
FEM_structure.LOAD_MATRIX=sparse(LOAD_MATRIX);
FEM_structure.Recovery_Matrix=sparse(Recovery_Matrix);
FEM_structure.U_static=U_static;
FEM_structure.FAN_Center_Recovery_vector=FAN_Center_Recovery_vector;
FEM_structure.U_FAN_static=U_FAN_static;
FEM_structure.COORD=COORD;
FEM_structure.ELEMENT=ELEMENT(:,2:9);
FEM_structure.NODE_SET=NODE_SET;
FEM_structure.Gamma=Gamma;
FEM_structure.intercoord=intercoord;
FEM_structure.FACES=FACES;
FEM_structure.ELint=ELint;
save FEM_s FEM_structure
clear
load FEM_s
print_all=false;
filter=true;
stress_plot=1;
FEM_structure.print_all=print_all;
Nthr=maxNumCompThreads(8);
%the point at the same x have the same 1 and 3 diplacement + linear
%% check if there are Lagrange multiplicator and eliminate orphan nodes
[FEM_structure] = No_lag(FEM_structure);
%% refinement
nfin=0;
[FEM_structure]=mesh_refinement_3D(FEM_structure,nfin);
%% Engine dofs eliminaion through fem shape function interpolation
%the point at the same x have the same 1 and 3 diplacement + linear
[FEM_structure]= engine_DOFs_elimination(FEM_structure);
%% Vectorized Stiffness matrix assembly
E0=210000;
[FEM_structure]= Stiffness_components3D(FEM_structure,E0);

%% INITIALIZE ITERATION
volfrac=0.15;
penal=4; ft=2;
x = volfrac*ones(size(FEM_structure.DESIGN_ZONE_ELEMENT_NODES,1),1);
xPhys = x(:);

% xPhys(x1>=3000&x4<=5000&z1<=min(z1)+dz)=0;
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
c       = 1000*eeem;
d       = eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 1000;
kkttol  = 1e-3;

%
%%%% The iterations start:
kktnorm = kkttol+1000;
outit = 0;
%% Preparing Filter
rmin=1;
[H,Hs]=preparing_Filter3D(FEM_structure,rmin);
%% FE-ANALYSIS
mgcg_selector=false;
[perfo,dperfo,GV,dGv,FEM_structure]=FEA3D(FEM_structure,xPhys,E0,penal,volfrac,mgcg_selector,nfin);
%% Test sensitivity with finite difference
test_sens=~true;
if test_sens
Ntest=10;  %number of tested variables
Range_test=randi(length(xPhys),Ntest,1);
ERR_p=zeros(Ntest,1);

dperfo0=dperfo;

x0=xPhys;
for sk=1:Ntest
k=Range_test(sk); 
pert=zeros(size(x0));
pert(k)=1e-4;
xp=x0+pert/2;
xm=x0-pert/2;
[perfop,~,~,~,FEM_structure]=FEA3D(FEM_structure,xp,E0,penal,volfrac,mgcg_selector,nfin);
[perfom,~,~,~,FEM_structure]=FEA3D(FEM_structure,xm,E0,penal,volfrac,mgcg_selector,nfin);
ERR_p(sk)=abs((perfop-perfom)/pert(k)-dperfo0(k))/(max(1e-4,abs((perfop-perfom)/pert(k))))*100;
end
figure
plot(ERR_p)
legend('Relative error on \Delta TSFC gradient (%)')
grid on
end
%% FILTERING/MODIFICATION OF SENSITIVITIES
if ft == 1
    dperfo(:) = H*(x(:).*dperfo.')./Hs./max(1e-3,x(:));
    %     dfandisp(:)=H*(x(:).*dfandisp.')./Hs./max(1e-3,x(:));
elseif ft == 2
    dperfo = H*(dperfo.'./Hs);
    %     dfandisp=H*(dfandisp.'./Hs);
    dGv(:) = H*(dGv(:)./Hs);
end
%% Gather results
f0val=1*perfo;
fval=GV;%
df0dx=1*reshape(dperfo,[],1);
dfdx=[dGv(:)'];%
innerit=0;
outvector1 = [outeriter innerit f0val fval'];
outvector2 = xval;
FACES=zeros(6*length(xval),4);
Fac_color=zeros(6*length(xval),1);
FID=[1 2 3 4;5 6 7 8;2 3 7 6;3 4 8 7;4 1 5 8;1 2 6 5];
for nn=1:length(xval)
    FACES(6*(nn-1)+(1:6),1)=FEM_structure.ELEMENT(nn,FID(:,1));
    FACES(6*(nn-1)+(1:6),2)=FEM_structure.ELEMENT(nn,FID(:,2));
    FACES(6*(nn-1)+(1:6),3)=FEM_structure.ELEMENT(nn,FID(:,3));
    FACES(6*(nn-1)+(1:6),4)=FEM_structure.ELEMENT(nn,FID(:,4));
    Fac_color(6*(nn-1)+(1:6))=(1-xval(nn))*ones(6,1);
end
%color filter
FACES_shown=FACES(Fac_color<=0.7,:);
Fac_color=Fac_color(Fac_color<=0.7);
h = figure(1); set(h,'Color',[1 1 1]);
clf
hold on; h24_patch = patch('Vertices',FEM_structure.COORD,'Faces',FACES_shown,'FaceVertexCData',Fac_color*[1 1 1],'FaceColor','flat'); axis equal; axis off;
quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
text(1000,0,0,'x')
text(0,1000,0,'y')
text(0,0,1000,'z')
view([27.6 18]);
drawnow;
hold off
fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',outeriter,perfo, ...
    GV,kktnorm);
%
figure(2)
hold on
scatter(outeriter,f0val,'fill','k')
figure(3)
% hold on
scatter(outeriter,GV,'fill','b')
folder=f0val;
error_old=100;
error=100;
%% START ITERATION
while kktnorm > kkttol & outit < maxoutit & error_old>0.5 & error>0.5
    outit   = outit+1;
    outeriter = outeriter+1;
    fold=f0val;
    %% MMA code optimization
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;
    if ft == 1
        xPhys=xval;
    elseif ft == 2
        xPhys = (H*xval(:))./Hs;
    end
    %     xPhys(x1>=2500&x4<=4500&z1==min(z1))=0;
    %% FE-ANALYSIS
    if rem(outit,15)==0
        Emin=FEM_structure.Emin0/outit;
    end
    mgcg_selector=false;
    [perfo,dperfo,GV,dGv,FEM_structure]=FEA3D(FEM_structure,xPhys,E0,penal,volfrac,mgcg_selector,nfin);
    
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dperfo(:) = H*(x(:).*dperfo.')./Hs./max(1e-3,x(:));
        %     dfandisp(:)=H*(x(:).*dfandisp.')./Hs./max(1e-3,x(:));
    elseif ft == 2
        dperfo = H*(dperfo.'./Hs);
        %     dfandisp=H*(dfandisp.'./Hs);
        dGv(:) = H*(dGv(:)./Hs);
    end
    %% Gather results
    f0val=1*perfo;
    fval=GV;%
    df0dx=1*reshape(dperfo,[],1);
    dfdx=[dGv(:)'];%
    innerit=0;
    outvector1 = [outeriter innerit f0val fval'];
    outvector2 = xval;
    Fac_color=zeros(6*length(xval),1);
    for nn=1:length(xval)
        Fac_color(6*(nn-1)+(1:6))=(1-xval(nn))*ones(6,1);
    end
    Optistruc_affich=0.7;
    Maff=1-Optistruc_affich;
    %color filter
    FACES_shown=FACES(Fac_color<=Maff,:);
    Fac_color=Fac_color(Fac_color<=Maff);
    h = figure(1); set(h,'Color',[1 1 1]);
    clf
    hold on; h24_patch = patch('Vertices',FEM_structure.COORD,'Faces',FACES_shown,'FaceVertexCData',Fac_color*[1 1 1],'FaceColor','flat'); axis equal; axis off;
    quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
    quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
    quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
    text(1000,0,0,'x')
    text(0,1000,0,'y')
    text(0,0,1000,'z')
    view([27.6 18]);
    drawnow;
    hold off
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',outeriter,perfo, ...
        GV,kktnorm);
    %
    figure(2)
    hold on
    scatter(outeriter,f0val,'fill','k')
    figure(3)
    % hold on
    scatter(outeriter,GV,'fill','b')
    %% %% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
    outvector1 = [outeriter innerit f0val fval(:)'];
    outvector2 = xval';
    error=abs((fold-f0val))/f0val*100;
    error_old=error;
end
% % Average dimension
% for nn=1:size(disp_DOF,1)
%     Grid_motion(new_DOFs_MAP(disp_DOF(nn),1),new_DOFs_MAP(disp_DOF(nn),3))=U(disp_DOF(nn));
% end
% MeanSize = mean(max(COORD)-min(COORD));
% % Maximum motion
% [MaxGridMotion,I] = max(abs(real(Grid_motion(:))));
% % MaxGridMotion = MaxGridMotion*sign(real(Grid_motion(I)));
% % New grid location
% %     NewGridPositionResidual = COORD(:,2:4) + 0.5*MeanSize*real(GridMotion_residual)/MaxGridMotionResidual;
% % Color displacement
% ColorOdsR = sqrt(sum(real(Grid_motion)'.^2))';
% ColorOdsR = ColorOdsR - min(ColorOdsR); ColorOdsR = ColorOdsR/max(ColorOdsR);
% % Plot
% 
% Final_coord=zeros(size(DESIGN_ZONE_COORD));
% for k=1:size(DESIGN_ZONE_COORD,1)
%     Final_coord(k,:)=DESIGN_ZONE_COORD(k,:)+0.5*MeanSize*[0,U(3*(k-1)+1),U(3*(k-1)+2),U(3*k)]/MaxGridMotion;
% end
% figure
% hold on; h25_patch = patch('Vertices',Final_coord(:,2:end),'Faces',FACES_shown(:,2:end),'CData',ColorOdsR,'FaceColor','interp','SpecularColorReflectance',0.3);
% colormap jet;
% quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
% quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
% quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
% text(1000,0,0,'x')
% text(0,1000,0,'y')
% text(0,0,1000,'z')
% view([27.6 18]);

