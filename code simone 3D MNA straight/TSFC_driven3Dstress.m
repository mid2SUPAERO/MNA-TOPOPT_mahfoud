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
print_all=true;
filter=true;
save_res=true;
FEM_structure.save_res=save_res;
FEM_structure.stress_plot=1;
FEM_structure.print_all=print_all;
% Nthr=maxNumCompThreads(8);
Engine_master=true;
plot_problem=false;
%the point at the same x have the same 1 and 3 diplacement + linear
%% check if there are Lagrange multiplicator and eliminate orphan nodes
[FEM_structure] = No_lag(FEM_structure);
%% refinement
nfin=2;
tic
[FEM_structure]=mesh_refinement_3D(FEM_structure,nfin);
toc
%% Engine dofs eliminaion through RBF interpolation
%the point at the same x have the same 1 and 3 diplacement + linear
tic
[FEM_structure]= engine_DOFs_eliminationpar(FEM_structure);
toc
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
FEM_structure.FACES=FACES;
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
tic
E0=210000;
[FEM_structure]= Stiffness_components3D(FEM_structure,E0);
toc
%% INITIALIZE ITERATION
volfrac=0.1;
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
kkttol  = 1e-2;

%
%%%% The iterations start:
kktnorm = kkttol+1000;
outit = 0;
%% Preparing Filter
rmin=1;
tic
[H,Hs]=preparing_Filter3Dpar(FEM_structure,rmin);
toc
%% FE-ANALYSIS
VMl=25;Penalty=4;
mgcg_selector=false;
Fac_color=zeros(6*length(xval),1);
for nn=1:length(xval)
    Fac_color(6*(nn-1)+(1:6))=(1-xval(nn))*ones(6,1);
end
Optistruc_affich=0.3;
Maff=1-Optistruc_affich;
%color filter
FACES_shown=FACES(Fac_color<=Maff,:);
Fac_color=Fac_color(Fac_color<=Maff);
toll_minres_U=1e-5;toll_minres_TSFC=1e-5;toll_minres_GKS=1e-5;
FEM_structure.FACES_shown=FACES_shown;
FEM_structure.outit=outit;
[perfo,dperfo,GV,dGv,GKSl,dGKSl,FEM_structure]=FEA3D(FEM_structure,xPhys,E0,penal,volfrac,mgcg_selector,nfin,Engine_master,VMl,Penalty,toll_minres_U,toll_minres_TSFC,toll_minres_GKS);
%% Test sensitivity with finite difference
test_sens=~true;
tollerance_U=1e-2.^[1:5];
nu=length(tollerance_U);
tollerance_gradient=1e-2.^[1:5];
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
        ERR_p_s=reshape(ERR_p(sk,:),nu,ng);
        ERR_G_s=reshape(ERR_p(sk,:),nu,ng);
        TU_s=reshape(TU(:),nu,ng);TG_s=reshape(TG(:),nu,ng);
        figure(20)
        hold on
        surf(TU_s,TG_s,ERR_p_s)
        xlabel('tolerance U')
        ylabel('tolerance on adjoint TSFC')
        zlabel('Relative error on \Delta TSFC% gradient (%)')
        figure(30)
        hold on
        surf(TU_s,TG_s,ERR_G_s)
        xlabel('tolerance U')
        ylabel('tolerance on adjoint GKS')
        zlabel('Relative error on GKs gradient (%)')
    end
    figure(40)
    ERR_p_s=reshape(max(ERR_p),nu,ng);
    ERR_G_s=reshape(max(ERR_G),nu,ng);
    hold on
    surf(TU_s,TG_s,ERR_p_s)
    xlabel('tolerance U')
    ylabel('tolerance on adjoint TSFC')
    zlabel('Relative error on \Delta TSFC% gradient (%)')
    figure(50)
    hold on
    surf(TU_s,TG_s,ERR_G_s)
    xlabel('tolerance U')
    ylabel('tolerance on adjoint GKS')
    zlabel('Relative error on GKs gradient (%)')
    %     plot(1:Ntest,ERR_p,'^-b',1:Ntest,ERR_G,'o-k')
    %     legend('Relative error on \Delta TSFC gradient (%)','Relative error on GKs gradient (%)')
    %     grid on
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
f0val=perfo;
fval=10*[GV;GKSl];%;1*(-0.38-FAN_Ax_disp)/0.38
df0dx=1*reshape(dperfo,[],1);
dfdx=10*[dGv(:)';dGKSl(:)'];%;-1*dfandisp(:)'/0.38
innerit=0;
outvector1 = [outeriter innerit f0val fval'];
outvector2 = xval;
if ~save_res
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
        print(['configuration_',num2str(outit,'%03d')],'-dpng')
    end
else
    coordo=FEM_structure.COORD;
    save(['configuration_',num2str(outit,'%03d')],'coordo','FACES_shown')
    clear coordo
end
error_old=100;
error=100;
fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f GKSl.:%7.3f kktn.:%7.3f och.:%7.3f\n',outeriter,perfo, ...
    Vol,GKSl,kktnorm,error);
%
if ~save_res
    figure(2)
    hold on
    plot(outeriter,perfo,'bo','MarkerFaceColor','b')
    plot(outeriter,Vol*100,'ro','MarkerFaceColor','r')
    plot(outeriter,(GKSl)*100,'ko','MarkerFaceColor','k')
    title(['Convergence volfrac = ',num2str(Vol*100),' %, \Delta TSFC =',num2str(perfo),'% ,G_{KS} \times 100 = ',num2str((GKSl)*100),' % , iter = ', num2str(outit)])
    grid on
    xlabel('iter')
    legend('\Delta TSFC %','Volume Fraction %','G_{KS} \times 100')
    if print_all
        print(['convergence_',num2str(outit,'%03d')],'-dpng')
    end
else
    save(['convergence_',num2str(outit,'%03d')],'outeriter','perfo','GKSl','Vol')
end
% figure(3)
% % hold on
% scatter(outeriter,GV,'fill','b')
folder=f0val;

%% START ITERATION
while  outit < maxoutit & ((error_old>0.05 | error>0.05) | GKSl>0.05|kktnorm > kkttol )
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
    if rem(outit,2000)==0
        FEM_structure.Emin=FEM_structure.Emin0/10^(fix(outit/20));
        penal=penal*1.05;
    end
    %% FE-ANALYSIS
    mgcg_selector=false;
    Fac_color=zeros(6*length(xval),1);
    for nn=1:length(xval)
        Fac_color(6*(nn-1)+(1:6))=(1-xval(nn))*ones(6,1);
    end
    Optistruc_affich=0.3;
    Maff=1-Optistruc_affich;
    %color filter
    FACES_shown=FACES(Fac_color<=Maff,:);
    Fac_color=Fac_color(Fac_color<=Maff);
    FEM_structure.outit=outit;
    FEM_structure.FACES_shown=FACES_shown;
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
    f0val=perfo;
    fval=10*[GV;GKSl];%;1*(-0.38-FAN_Ax_disp)/0.38
    df0dx=1*reshape(dperfo,[],1);
    dfdx=10*[dGv(:)';dGKSl(:)'];%;-1*dfandisp(:)'/0.38
    
    if ~save_res
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
            print(['configuration_',num2str(outit,'%03d')],'-dpng')
        end
    else
        save(['configuration_',num2str(outit,'%03d')],'FACES_shown')
    end
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f GKSl.:%7.3f kktn.:%7.3f och.:%7.3f\n',outeriter,perfo, ...
        Vol,GKSl,kktnorm,error);
    %
    if ~save_res
        figure(2)
        hold on
        plot(outeriter,perfo,'bo','MarkerFaceColor','b')
        plot(outeriter,Vol*100,'ro','MarkerFaceColor','r')
        plot(outeriter,(GKSl)*100,'ko','MarkerFaceColor','k')
        title(['Convergence volfrac = ',num2str(Vol*100),' %, \Delta TSFC =',num2str(perfo),'% ,G_{KS} \times 100 = ',num2str((GKSl)*100),' % , iter = ', num2str(outit)])
        grid on
        xlabel('iter')
        legend('\Delta TSFC %','Volume Fraction %','G_{KS} \times 100')
        if print_all
            print(['convergence_',num2str(outit,'%03d')],'-dpng')
        end
    else
       save(['convergence_',num2str(outit,'%03d')],'outeriter','perfo','GKSl','Vol')
    end
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

