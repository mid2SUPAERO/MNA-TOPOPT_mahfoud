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
print_all=~true;
filter=true;
stress_plot=1;
FEM_structure.print_all=print_all;
Nthr=maxNumCompThreads(8);
Engine_master=true;
plot_problem=~false;
%the point at the same x have the same 1 and 3 diplacement + linear
%% check if there are Lagrange multiplicator and eliminate orphan nodes
[FEM_structure] = No_lag(FEM_structure);
%% refinement
nfin=1;
[FEM_structure]=mesh_refinement_3D(FEM_structure,nfin);
%% Engine dofs eliminaion through fem shape function interpolation
%the point at the same x have the same 1 and 3 diplacement + linear
[FEM_structure]= engine_DOFs_elimination(FEM_structure);
%% plot the engine and DZ mesh together
El_numb=size(FEM_structure.ELEMENT,1);
FACES=zeros(6*El_numb,4);
FID=[1 2 3 4;5 6 7 8;2 3 7 6;3 4 8 7;4 1 5 8;1 2 6 5];
for nn=1:El_numb
    FACES(6*(nn-1)+(1:6),1)=FEM_structure.ELEMENT(nn,FID(:,1));
    FACES(6*(nn-1)+(1:6),2)=FEM_structure.ELEMENT(nn,FID(:,2));
    FACES(6*(nn-1)+(1:6),3)=FEM_structure.ELEMENT(nn,FID(:,3));
    FACES(6*(nn-1)+(1:6),4)=FEM_structure.ELEMENT(nn,FID(:,4));
end
if plot_problem
    %%show only external faces for clarity
    % Id=zeros(6*El_numb,1);
    % for nn=1:6*El_numb
    %     Id(nn)=sum(sum(ismember(FACES,FACES(nn,:)),2)==4)==1;
    % end
    FACE_intE=FEM_structure.ELint;
    for k=1:size(FEM_structure.ELint,1)
        for l=1:size(FEM_structure.ELint,2)-1
            FACE_intE(k,l+1)=find(FEM_structure.intercoord(:,1)==FEM_structure.ELint(k,l+1));
        end
    end
    h = figure(3);clf
    set(h,'Color',[1 1 1]);
    axes1 = axes('Parent',h);
    hold on; h24_patch = patch('Parent',axes1,'Vertices',FEM_structure.COORD,'Faces',FACES,'FaceVertexCData',[.5 .1 .1],'FaceColor','flat'); axis equal; axis off;
    hold on; h25_patch = patch('Parent',axes1,'Vertices',FEM_structure.intercoord(:,2:4),'Faces',FACE_intE(:,2:5),'FaceVertexCData',[.6 .5 .3],'FaceColor','flat'); axis equal; axis off;
    s1=quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k');
    s1=quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k');
    s1=quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k');
    text(1000,0,0,'x')
    text(0,1000,0,'y')
    text(0,0,1000,'z')
    view([27.6 18]);
    legend([h24_patch,h25_patch],'Design Zone FEM','Retained Engine Interface')
    % Set the remaining axes properties
    set(axes1,'CameraViewAngle',3.84174456530832,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatio',[434 342.3 684.6]);
    drawnow;
end

%% Vectorized Stiffness matrix assembly
E0=210000;
[FEM_structure]= Stiffness_components3D(FEM_structure,E0);

%% INITIALIZE ITERATION
volfrac=0.2;
admissible_DTSFC=2;
penal=5; ft=2;
% x = volfrac*ones(size(FEM_structure.DESIGN_ZONE_ELEMENT_NODES,1),1);
x = ones(size(FEM_structure.DESIGN_ZONE_ELEMENT_NODES,1),1);

xPhys = x(:);

% xPhys(x1>=3000&x4<=5000&z1<=min(z1)+dz)=0;
loop = 0;
change = 1;
m = 2;
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
VMl=25;Penalty=4;
mgcg_selector=false;
toll_minres_U=1e-5;toll_minres_TSFC=1e-5;toll_minres_GKS=1e-5;
[perfo,dperfo,GV,dGv,GKSl,dGKSl,FEM_structure]=FEA3D(FEM_structure,xPhys,E0,penal,volfrac,mgcg_selector,nfin,Engine_master,VMl,Penalty,toll_minres_U,toll_minres_TSFC,toll_minres_GKS);
%% Test sensitivity with finite difference
test_sens=~true;
tollerance_U=1e-1.^[5];
nu=length(tollerance_U);
tollerance_gradient=1e-1.^[5];
ng=length(tollerance_gradient);
[TU,TG]=meshgrid(tollerance_U,tollerance_gradient);TU=TU(:);TG=TG(:); ntol=length(TU);

if test_sens
    Ntest=10;  %number of tested variables
    Range_test=randi(length(xPhys),Ntest,1);
    ERR_p=zeros(Ntest,ntol);
    ERR_G=zeros(Ntest,ntol);
    dperfo0=dperfo;
    dGKSl0=dGKSl;
    x0=xPhys;
    for sk=1:Ntest
        k=Range_test(sk);
        pert=zeros(size(x0));
        pert(k)=1e-4;
        xp=x0+pert/2;
        xm=x0-pert/2;
        for ki_tol=1:ntol
            FEM_structure.Sol=zeros(size(FEM_structure.Sol));
            FEM_structure.psi_tild_sol=zeros(size(FEM_structure.psi_tild_sol));
            [perfop,~,~,~,GKSlp,~,FEM_structure]=FEA3D(FEM_structure,xp,E0,penal,volfrac,mgcg_selector,nfin,Engine_master,VMl,Penalty,TU(ki_tol),TG(ki_tol),TG(ki_tol));
            FEM_structure.Sol=zeros(size(FEM_structure.Sol));
            FEM_structure.psi_tild_sol=zeros(size(FEM_structure.psi_tild_sol));
            [perfom,~,~,~,GKSlm,~,FEM_structure]=FEA3D(FEM_structure,xm,E0,penal,volfrac,mgcg_selector,nfin,Engine_master,VMl,Penalty,TU(ki_tol),TG(ki_tol),TG(ki_tol));
            ERR_p(sk,ki_tol)=abs((perfop-perfom)/pert(k)-dperfo0(k))/(max(1e-4,abs((perfop-perfom)/pert(k))))*100;
            ERR_G(sk,ki_tol)=abs((GKSlp-GKSlm)/pert(k)-dGKSl0(k))/(max(1e-4,abs((GKSlp-GKSlm)/pert(k))))*100;
        end
    end
    figure
    plot(1:Ntest,ERR_p,'^-b',1:Ntest,ERR_G,'o-k')
    legend('Relative error on \Delta TSFC gradient (%)','Relative error on GKs gradient (%)')
    grid on
    %% convergence test
    [perfo_r,dperfo_r,~,~,GKSl_r,dGKSl_r,FEM_structure]=FEA3D(FEM_structure,xPhys,E0,penal,volfrac,mgcg_selector,nfin,Engine_master,VMl,Penalty,0,0,0);
    
    tollerance_U=1e-1;tollerance_gradient=tollerance_U;
    for k=1:10
        FEM_structure.Sol=zeros(size(FEM_structure.Sol));
        FEM_structure.psi_tild_sol=zeros(size(FEM_structure.psi_tild_sol));
        [perfo_e(k),dperfo_e(k,:),~,~,GKSl_e(k),dGKSl_e(k,:),FEM_structure]=FEA3D(FEM_structure,xPhys,E0,penal,volfrac,mgcg_selector,nfin,Engine_master,VMl,Penalty,tollerance_U^k,tollerance_gradient^k,tollerance_gradient^k);
    end
    Ep=(perfo_e-perfo_r)/perfo_r*100;
    Egp=sqrt(sum((dperfo_e-repmat(dperfo_r,10,1)).^2,2))/norm(dperfo_r)*100;
    Eg=(GKSl_e-GKSl_r)/GKSl_r*100;
    Egg=sqrt(sum((dGKSl_e-repmat(dGKSl_r,10,1)).^2,2))/norm(dGKSl_r)*100;
    figure
    semilogx(tollerance_U.^(1:10),Ep,'-ok','LineWidth',2,'MarkerFaceColor','k')
    hold on
    semilogx(tollerance_U.^(1:10),Egp,'-or','LineWidth',2,'MarkerFaceColor','r')
    semilogx(tollerance_U.^(1:10),Eg,'-om','LineWidth',2,'MarkerFaceColor','m')
    semilogx(tollerance_U.^(1:10),Egg,'-ob','LineWidth',2,'MarkerFaceColor','b')
    grid on
    legend('\Delta TSFC % relative error [%]','\Delta TSFC % gradient relative error [%]','G_{KS} relative error [%]','G_{KS} gradient relative error [%]')
    xlabel('minres residual tolerance')
    
end
%% FILTERING/MODIFICATION OF SENSITIVITIES
if ft == 1
    dperfo(:) = H*(x(:).*dperfo(:))./Hs./max(1e-3,x(:));
    dGKSl(:)= H*(x(:).*dGKSl(:))./Hs./max(1e-3,x(:));
    %     dfandisp(:)=H*(x(:).*dfandisp.')./Hs./max(1e-3,x(:));
elseif ft == 2
    dperfo = H*(dperfo(:)./Hs);
    dGKSl(:)= H*(dGKSl(:)./Hs);
    %     dfandisp=H*(dfandisp.'./Hs);
    dGv(:) = H*(dGv(:)./Hs);
end
%% Gather results
Gp=(perfo-admissible_DTSFC)/admissible_DTSFC;
dGp=dperfo/admissible_DTSFC;
Vol=(1+GV)*volfrac;
dVol=dGv*volfrac;
f0val=Vol;
fval=10*[Gp;GKSl];%;1*(-0.38-FAN_Ax_disp)/0.38
df0dx=1*reshape(dVol,[],1);
dfdx=10*[dGp(:)';dGKSl(:)'];%;-1*dfandisp(:)'/0.38
innerit=0;
outvector1 = [outeriter innerit f0val fval'];
outvector2 = xval;
Fac_color=zeros(6*length(xval),1);
for nn=1:length(xval)
    Fac_color(6*(nn-1)+(1:6))=(1-xval(nn))*ones(6,1);
end
%color filter
FACES_shown=FACES(Fac_color<=0.7,:);
Fac_color=Fac_color(Fac_color<=0.7);
h = figure(1);clf
set(h,'Color',[1 1 1]);
axes1 = axes('Parent',h);
hold on; h24_patch = patch('Parent',axes1,'Vertices',FEM_structure.COORD,'Faces',FACES_shown,'FaceVertexCData',[.5 .1 .1],'FaceColor','flat'); axis equal; axis off;
quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
text(1000,0,0,'x')
text(0,1000,0,'y')
text(0,0,1000,'z')
view([27.6 18]);
% Set the remaining axes properties
set(axes1,'CameraViewAngle',3.84174456530832,'DataAspectRatio',[1 1 1],...
    'PlotBoxAspectRatio',[434 342.3 684.6]);
drawnow;
hold off
if print_all
    print(['configuration_',num2str(outit,'%03d')],'-djpeg')
end
error_old=100;
error=100;
fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f GKSl.:%7.3f kktn.:%7.3f och.:%7.3f\n',outeriter,perfo, ...
    Vol,GKSl,kktnorm,error);
%
figure(2)
hold on
plot(outeriter,perfo,'bo','MarkerFaceColor','b')
plot(outeriter,Vol*100,'ro','MarkerFaceColor','r')
plot(outeriter,(1+GKSl)*VMl,'ko','MarkerFaceColor','k')
title(['Convergence volfrac = ',num2str(Vol*100),', \Delta TSFC% =',num2str(perfo),',\sigma_{max} \approx ',num2str((1+GKSl)*VMl),' MPa , iter = ', num2str(outit)])
grid on
xlabel('iter')
legend('\Delta TSFC %','Volume Fraction %','\sigma_{max}')
if print_all
    print(['convergence_',num2str(outit,'%03d')],'-djpeg')
end
% figure(3)
% % hold on
% scatter(outeriter,GV,'fill','b')
folder=f0val;
x1=FEM_structure.x1;
y1=FEM_structure.y1;
z1=FEM_structure.z1;
x2=FEM_structure.x2;
y2=FEM_structure.y2;
z2=FEM_structure.z2;
x3=FEM_structure.x3;
y3=FEM_structure.y3;
z3=FEM_structure.z3;
x4=FEM_structure.x4;
y4=FEM_structure.y4;
z4=FEM_structure.z4;
x5=FEM_structure.x5;
y5=FEM_structure.y5;
z5=FEM_structure.z5;
x6=FEM_structure.x6;
y6=FEM_structure.y6;
z6=FEM_structure.z6;
x7=FEM_structure.x7;
y7=FEM_structure.y7;
z7=FEM_structure.z7;
x8=FEM_structure.x8;
y8=FEM_structure.y8;
z8=FEM_structure.z8;
center_node_coordinate=[(x1+x2+x3+x4+x5+x6+x7+x8)/8,(y1+y2+y3+y4+y5+y6+y7+y8)/8,(z1+z2+z3+z4+z5+z6+z7+z8)/8];
[center_node_coordinate,Ic]=sortrows(center_node_coordinate);
Xg=reshape(center_node_coordinate(:,1),8*(nfin+1),25*(nfin+1),38*(nfin+1));
Yg=reshape(center_node_coordinate(:,2),8*(nfin+1),25*(nfin+1),38*(nfin+1));
Zg=reshape(center_node_coordinate(:,3),8*(nfin+1),25*(nfin+1),38*(nfin+1));
%% START ITERATION
while kktnorm > kkttol & outit < maxoutit & ((error_old>0.005 | error>0.005) | (Gp>0.005)| GKSl>0.005)
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
    %% penalty correction
    if rem(outit,200)==0
        Emin=FEM_structure.Emin0/outit;
        penal=penal*1.05;
    end
    %% FE-ANALYSIS
    VMl=25;Penalty=4;
    mgcg_selector=false;
    [perfo,dperfo,GV,dGv,GKSl,dGKSl,FEM_structure]=FEA3D(FEM_structure,xPhys,E0,penal,volfrac,mgcg_selector,nfin,Engine_master,VMl,Penalty,toll_minres_U,toll_minres_TSFC,toll_minres_GKS);
    
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dperfo(:) = H*(x(:).*dperfo(:))./Hs./max(1e-3,x(:));
        dGKSl(:)= H*(x(:).*dGKSl(:))./Hs./max(1e-3,x(:));
        %     dfandisp(:)=H*(x(:).*dfandisp.')./Hs./max(1e-3,x(:));
    elseif ft == 2
        dperfo = H*(dperfo(:)./Hs);
        dGKSl(:)= H*(dGKSl(:)./Hs);
        %     dfandisp=H*(dfandisp.'./Hs);
        dGv(:) = H*(dGv(:)./Hs);
    end
    %% Gather results
    Gp=(perfo-admissible_DTSFC)/admissible_DTSFC;
    dGp=dperfo/admissible_DTSFC;
    Vol=(1+GV)*volfrac;
    dVol=dGv*volfrac;
    f0val=Vol;
    fval=10*[Gp;GKSl];%;1*(-0.38-FAN_Ax_disp)/0.38
    df0dx=1*reshape(dVol,[],1);
    dfdx=10*[dGp(:)';dGKSl(:)'];%;-1*dfandisp(:)'/0.38
    Fac_color=zeros(6*length(xval),1);
    for nn=1:length(xval)
        Fac_color(6*(nn-1)+(1:6))=(1-xval(nn))*ones(6,1);
    end
    Optistruc_affich=0.3;
    Maff=1-Optistruc_affich;
    %color filter
    FACES_shown=FACES(Fac_color<=Maff,:);
    Fac_color=Fac_color(Fac_color<=Maff);
    h = figure(1);clf
    set(h,'Color',[1 1 1]);
    axes1 = axes('Parent',h);
    hold on; h24_patch = patch('Parent',axes1,'Vertices',FEM_structure.COORD,'Faces',FACES_shown,'FaceVertexCData',[.5 .1 .1],'FaceColor','flat'); axis equal; axis off;
    quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
    quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
    quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
    text(1000,0,0,'x')
    text(0,1000,0,'y')
    text(0,0,1000,'z')
    view([27.6 18]);
    % Set the remaining axes properties
    set(axes1,'CameraViewAngle',3.84174456530832,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatio',[434 342.3 684.6]);
    drawnow;
    hold off
    if print_all
        print(['configuration_',num2str(outit,'%03d')],'-djpeg')
    end
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f GKSl.:%7.3f kktn.:%7.3f och.:%7.3f\n',outeriter,perfo, ...
        Vol,GKSl,kktnorm,error);
    %
    figure(2)
    hold on
    plot(outeriter,perfo,'bo','MarkerFaceColor','b')
    plot(outeriter,Vol*100,'ro','MarkerFaceColor','r')
    plot(outeriter,(1+GKSl)*VMl,'ko','MarkerFaceColor','k')
    title(['Convergence volfrac = ',num2str(Vol*100),', \Delta TSFC% =',num2str(perfo),',\sigma_{max} \approx ',num2str((1+GKSl)*VMl),' MPa , iter = ', num2str(outit)])
    grid on
    xlabel('iter')
    legend('\Delta TSFC %','Volume Fraction %','\sigma_{max}')
    if print_all
        print(['convergence_',num2str(outit,'%03d')],'-djpeg')
    end
    figure(8)
    clf
    hold on
    Xs=reshape(xPhys(Ic),8*(nfin+1),25*(nfin+1),38*(nfin+1));
    p = patch(isosurface(Xg,Yg,Zg,Xs,0.5));
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    camlight
    lighting gouraud
    quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
    quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
    quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
    text(1000,0,0,'x')
    text(0,1000,0,'y')
    text(0,0,1000,'z')
    view([27.6 18]);
    axis equal
    axis off
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

