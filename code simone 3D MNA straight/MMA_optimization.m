function [C,perfo,Vol,x,kktnorm]=MMA_optimization(FEM_structure,formulation)
%% Path generation
switch  formulation.type
    case 1
        objective_string='C';
        constraint_strinc=['V_',num2str(round(formulation.volfrac*100)),'_T_',num2str(round(formulation.admissible_DTSFC*1000))];
    case 2
        objective_string='V';
        constraint_strinc=['C_',num2str(round(formulation.admissible_compliance)),'_T_',num2str(round(formulation.admissible_DTSFC*1000))];
    case 3
        objective_string='T';
        constraint_strinc=['V_',num2str(round(formulation.volfrac*100)),'_C_',num2str(round(formulation.admissible_compliance))];
end
Path=['min_',objective_string,'st_',constraint_strinc,'/'];
mkdir(['min_',objective_string,'st_',constraint_strinc])
%% FE-ANALYSIS
FEM_structure.VMl=10;FEM_structure.Aggregation_constant=4;
disp(['Maximum stress allowed ',num2str(FEM_structure.VMl*(1+log(8*FEM_structure.n)/FEM_structure.Aggregation_constant))])
mgcg_selector=false;
Fac_color=zeros(6*length(formulation.xPhys),1); % TO DO: update?
for nn=1:length(formulation.xPhys) % TO DO: update?
    Fac_color(6*(nn-1)+(1:6))=(1-formulation.xPhys(nn))*ones(6,1); % TO DO: update?
end
Optistruc_affich=0.3;
Maff=1-Optistruc_affich;
%color filter
FACES_shown=FEM_structure.FACES(Fac_color<=Maff,:);
Fac_color=Fac_color(Fac_color<=Maff);
FEM_structure.toll_minres_U=1e-5;FEM_structure.toll_minres_TSFC=1e-5;FEM_structure.toll_minres_GKS=1e-5;
FEM_structure.FACES_shown=FACES_shown;
FEM_structure.FEM_structure.outit=FEM_structure.outit;
[C,dC,perfo,dperfo,GV,dGv,GKSl,dGKSl,FEM_structure]=FEA3D_compliance(FEM_structure,formulation.xPhys,formulation.den_grad,FEM_structure.E0,FEM_structure.penal,formulation.volfrac,mgcg_selector,FEM_structure.nfin,FEM_structure.Engine_master,FEM_structure.VMl,FEM_structure.Aggregation_constant,FEM_structure.toll_minres_U,FEM_structure.toll_minres_TSFC,FEM_structure.toll_minres_GKS);
%% Test sensitivity with finite difference
test_sens=~true ;
tollerance_U=1e-1.^[5] ;
nu=length(tollerance_U) ;
tollerance_gradient=1e-1.^[5] ;
ng=length(tollerance_gradient) ;
[TU,TG]=meshgrid(tollerance_U,tollerance_gradient);TU=TU(:);TG=TG(:); ntol=length(TU);

if test_sens 
    Ntest=10;  %number of tested variables
    Range_test=randi(length(formulation.x),Ntest,1);
    ERR_p=zeros(Ntest,ntol);
    ERR_G=zeros(Ntest,ntol);
    dperfo0=dperfo;
    dGKSl0=dGKSl;
    x0=formulation.x;
    error = zeros(nn,length(x0));
    for sk=1:Ntest
        k=Range_test(sk);
        pert=zeros(size(x0));
        pert(k)=1e-4;
        xp=x0+pert/2;
        xm=x0-pert/2;
        [xp_d,den_gradp]=mna_density(FEM_structure,xp) ;
        [xm_d,den_gradm]=mna_density(FEM_structure,xm) ;
        error(:,k)=abs((xp_d-xm_d)/pert(k)-(den_gradp(:,k)+den_gradm(:,k))/2)./max(abs((den_gradp(:,k)+den_gradm(:,k))/2),1e-2) ;
        for ki_tol=1:ntol
            FEM_structure.Sol=zeros(size(FEM_structure.Sol));
            FEM_structure.psi_tild_sol=zeros(size(FEM_structure.psi_tild_sol));
            [~,~,perfop,~,~,~,GKSlp,~,FEM_structure]=FEA3D_compliance(FEM_structure,xp_d,den_gradp,FEM_structure.E0,FEM_structure.penal,formulation.volfrac,mgcg_selector,FEM_structure.nfin,FEM_structure.Engine_master,FEM_structure.VMl,FEM_structure.Aggregation_constant,TU(ki_tol),TG(ki_tol),TG(ki_tol));
            FEM_structure.Sol=zeros(size(FEM_structure.Sol));
            FEM_structure.psi_tild_sol=zeros(size(FEM_structure.psi_tild_sol));
            [~,~,perfom,~,~,~,GKSlm,~,FEM_structure]=FEA3D_compliance(FEM_structure,xm_d,den_gradm,FEM_structure.E0,FEM_structure.penal,formulation.volfrac,mgcg_selector,FEM_structure.nfin,FEM_structure.Engine_master,FEM_structure.VMl,FEM_structure.Aggregation_constant,TU(ki_tol),TG(ki_tol),TG(ki_tol));
            ERR_p(sk,ki_tol)=abs((perfop-perfom)/pert(k)-dperfo0(k))/(max(1e-4,abs((perfop-perfom)/pert(k))))*100;
            ERR_G(sk,ki_tol)=abs((GKSlp-GKSlm)/pert(k)-dGKSl0(k))/(max(1e-4,abs((GKSlp-GKSlm)/pert(k))))*100;
        end
    end
    figure
    plot(1:Ntest,ERR_p,'^-b',1:Ntest,ERR_G,'o-k')
    legend('Relative error on \Delta TSFC gradient (%)','Relative error on GKs gradient (%)')
    grid on
    %     %% convergence test
    %     [perfo_r,dperfo_r,~,~,GKSl_r,dGKSl_r,FEM_structure]=FEA3D(FEM_structure,formulation.xPhys,FEM_structure.E0,FEM_structure.penal,formulation.volfrac,mgcg_selector,FEM_structure.nfin,FEM_structure.Engine_master,FEM_structure.VMl,FEM_structure.Aggregation_constant,0,0,0);
    %     % FEA3D doesn't exist !!!!
    %     tollerance_U=1e-1;tollerance_gradient=tollerance_U;
    %     for k=1:10
    %         FEM_structure.Sol=zeros(size(FEM_structure.Sol));
    %         FEM_structure.psi_tild_sol=zeros(size(FEM_structure.psi_tild_sol));
    %         [perfo_e(k),dperfo_e(k,:),~,~,GKSl_e(k),dGKSl_e(k,:),FEM_structure]=FEA3D(FEM_structure,formulation.xPhys,FEM_structure.E0,FEM_structure.penal,formulation.volfrac,mgcg_selector,FEM_structure.nfin,FEM_structure.Engine_master,FEM_structure.VMl,FEM_structure.Aggregation_constant,tollerance_U^k,tollerance_gradient^k,tollerance_gradient^k);
    %     end
    %     Ep=(perfo_e-perfo_r)/perfo_r*100;
    %     Egp=sqrt(sum((dperfo_e-repmat(dperfo_r,10,1)).^2,2))/norm(dperfo_r)*100 ;
    %     Eg=(GKSl_e-GKSl_r)/GKSl_r*100 ;
    %     Egg=sqrt(sum((dGKSl_e-repmat(dGKSl_r,10,1)).^2,2))/norm(dGKSl_r)*100 ;
    %     figure
    %     semilogx(tollerance_U.^(1:10),Ep,'-ok','LineWidth',2,'MarkerFaceColor','k')
    %     hold on
    %     semilogx(tollerance_U.^(1:10),Egp,'-or','LineWidth',2,'MarkerFaceColor','r')
    %     semilogx(tollerance_U.^(1:10),Eg,'-om','LineWidth',2,'MarkerFaceColor','m')
    %     semilogx(tollerance_U.^(1:10),Egg,'-ob','LineWidth',2,'MarkerFaceColor','b')
    %     grid on
    %     legend('\Delta TSFC % relative error [%]','\Delta TSFC % gradient relative error [%]','G_{KS} relative error [%]','G_{KS} gradient relative error [%]')
    %     xlabel('minres residual tolerance')
    
end
% %% FILTERING/MODIFICATION OF SENSITIVITIES
% if FEM_structure.ft == 1 % TO DO: can be removed §§§
%     dperfo(:) = FEM_structure.H*(formulation.x(:).*dperfo(:))./FEM_structure.Hs./max(1e-3,formulation.x(:));
%     dGKSl(:)= FEM_structure.H*(formulation.x(:).*dGKSl(:))./FEM_structure.Hs./max(1e-3,formulation.x(:));
%     dC(:)= FEM_structure.H*(formulation.x(:).*dC(:))./FEM_structure.Hs./max(1e-3,formulation.x(:));
%     %     dfandisp(:)=H*(x(:).*dfandisp.')./Hs./max(1e-3,x(:));
% elseif FEM_structure.ft == 2
%     dperfo = FEM_structure.H*(dperfo(:)./FEM_structure.Hs);
%     dGKSl(:)= FEM_structure.H*(dGKSl(:)./FEM_structure.Hs);
%     %     dfandisp=H*(dfandisp.'./Hs);
%     dGv(:) = FEM_structure.H*(dGv(:)./FEM_structure.Hs);
%     dC(:) = FEM_structure.H*(dC(:)./FEM_structure.Hs);
% end
%% Gather results
Gc=(C-formulation.admissible_compliance)/formulation.admissible_compliance*100;
Gp=(perfo-formulation.admissible_DTSFC)/formulation.admissible_DTSFC*100;
dGp=dperfo/formulation.admissible_DTSFC*100;
dGc=dC/formulation.admissible_compliance*100;
Vol=(1+GV)*formulation.volfrac;
dVol=dGv*formulation.volfrac;
switch formulation.type
    case 1
        f0val=C;
        fval=100*[GV*100;GKSl;Gp];%;1*(-0.38-FAN_Ax_disp)/0.38
        df0dx=1*reshape(dC,[],1);
        dfdx=100*[dGv(:)'*100;dGKSl(:)';dGp(:)'];%;-1*dfandisp(:)'/0.38
    case 2
        f0val=Vol*100;
        fval=10*[Gc;GKSl;Gp];%;1*(-0.38-FAN_Ax_disp)/0.38
        df0dx=100*reshape(dVol,[],1);
        dfdx=10*[dGc(:)';dGKSl(:)';dGp(:)'];%;-1*dfandisp(:)'/0.38
    case 3
        f0val=perfo*100;
        fval=10*[Gc;GKSl;GV];%;1*(-0.38-FAN_Ax_disp)/0.38
        df0dx=100*reshape(dperfo,[],1);
        dfdx=10*[dGc(:)';dGKSl(:)';dGv(:)'];%;-1*dfandisp(:)'/0.38
end
innerit=0;
outvector1 = [FEM_structure.outeriter innerit f0val fval'];
outvector2 = formulation.xPhys; % TO DO: update?
Fac_color=zeros(6*length(formulation.xPhys),1); % TO DO: update?
for nn=1:length(formulation.xPhys) % TO DO: update?
    Fac_color(6*(nn-1)+(1:6))=(1-formulation.xPhys(nn))*ones(6,1); % TO DO: update?
end
%color filter
FACES_shown1=FEM_structure.FACES(Fac_color<=0.7,:);
Fac_color1=Fac_color(Fac_color<=0.7);
FACES_shown2=FEM_structure.FACES(Fac_color<=0.1,:);
Fac_color2=Fac_color(Fac_color<=0.1);
% h = figure(1);clf
% set(h,'Color',[1 1 1]);
% axes1 = axes('Parent',h);
% hold on; h24_patch = patch('Parent',axes1,'Vertices',FEM_structure.COORD,'Faces',FACES_shown1,'FaceVertexCData',Fac_color1,'FaceColor','flat'); axis equal; axis off;
% quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
% quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
% quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
% text(1000,0,0,'x')
% text(0,1000,0,'y')
% text(0,0,1000,'z')
% view([27.6 18]);
% % Set the remaining axes properties
% set(axes1,'CameraViewAngle',3.84174456530832,'DataAspectRatio',[1 1 1],...
%     'PlotBoxAspectRatio',[434 342.3 684.6]);
% drawnow;

% hold off
% if print_all
%     print(['configuration_',num2str(outit,'%03d')],'-djpeg')
% end
% h = figure(6);clf
% set(h,'Color',[1 1 1]);
% axes1 = axes('Parent',h);
% hold on; h24_patch = patch('Parent',axes1,'Vertices',FEM_structure.COORD,'Faces',FACES_shown2,'FaceVertexCData',Fac_color2,'FaceColor','flat'); axis equal; axis off;
% quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
% quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
% quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
% text(1000,0,0,'x')
% text(0,1000,0,'y')
% text(0,0,1000,'z')
% view([27.6 18]);
% Set the remaining axes properties
% set(axes1,'CameraViewAngle',3.84174456530832,'DataAspectRatio',[1 1 1],...
%     'PlotBoxAspectRatio',[434 342.3 684.6]);
% drawnow;
error_old=100;
error=100;
fprintf(' It.:%5i C:%11.4f T.:%11.4f Vol.:%7.3f GKSl.:%7.3f kktn.:%7.3f och.:%7.3f\n',FEM_structure.outit,C,perfo ,...
    Vol,GKSl,FEM_structure.kktnorm,error);
%
% figure(2)
% hold on
% plot(outeriter,perfo,'bo','MarkerFaceColor','b')
% plot(outeriter,Vol*100,'ro','MarkerFaceColor','r')
% plot(outeriter,(1+GKSl)*VMl,'ko','MarkerFaceColor','k')
% title(['Convergence volfrac = ',num2str(Vol*100),', \Delta TSFC% =',num2str(C),',\sigma_{max} \approx ',num2str((1+GKSl)*VMl),' MPa , iter = ', num2str(outit)])
% grid on
% xlabel('iter')
% legend('\Delta TSFC %','Volume Fraction %','\sigma_{max}')
% if print_all
%     print(['convergence_',num2str(outit,'%03d')],'-djpeg')
% end
outit=FEM_structure.outit;
save([Path,'convergence_',num2str(FEM_structure.outit,'%03d')],'outit','perfo','GKSl','Vol','C')
coordo=FEM_structure.COORD;
save(['configuration_',num2str(FEM_structure.outit,'%03d')],'coordo','FACES_shown')
clear coordo
% figure(3)
% % hold on
% scatter(outeriter,GV,'fill','b')
% folder=f0val;
% x1=FEM_structure.x1;
% y1=FEM_structure.y1;
% z1=FEM_structure.z1;
% x2=FEM_structure.x2;
% y2=FEM_structure.y2;
% z2=FEM_structure.z2;
% x3=FEM_structure.x3;
% y3=FEM_structure.y3;
% z3=FEM_structure.z3;
% x4=FEM_structure.x4;
% y4=FEM_structure.y4;
% z4=FEM_structure.z4;
% x5=FEM_structure.x5;
% y5=FEM_structure.y5;
% z5=FEM_structure.z5;
% x6=FEM_structure.x6;
% y6=FEM_structure.y6;
% z6=FEM_structure.z6;
% x7=FEM_structure.x7;
% y7=FEM_structure.y7;
% z7=FEM_structure.z7;
% x8=FEM_structure.x8;
% y8=FEM_structure.y8;
% z8=FEM_structure.z8;
% center_node_coordinate=[(x1+x2+x3+x4+x5+x6+x7+x8)/8,(y1+y2+y3+y4+y5+y6+y7+y8)/8,(z1+z2+z3+z4+z5+z6+z7+z8)/8];
% [center_node_coordinate,Ic]=sortrows(center_node_coordinate);
% Xg=reshape(center_node_coordinate(:,1),8*(nfin+1),25*(nfin+1),38*(nfin+1));
% Yg=reshape(center_node_coordinate(:,2),8*(nfin+1),25*(nfin+1),38*(nfin+1));
% Zg=reshape(center_node_coordinate(:,3),8*(nfin+1),25*(nfin+1),38*(nfin+1));
%% START ITERATION
change=10;
f3 = figure() ;
while FEM_structure.kktnorm > FEM_structure.kkttol && FEM_structure.outit < FEM_structure.maxoutit
    FEM_structure.outit   = FEM_structure.outit+1;
    FEM_structure.outeriter = FEM_structure.outeriter+1;
    fold=f0val;
%     mna3D_plot(formulation,FEM_structure) ;
    %% MMA code optimization
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,S,FEM_structure.low,FEM_structure.upp] = ...  % TO DO: update §§§
        mmasub(FEM_structure.m,FEM_structure.n_var,FEM_structure.outeriter,formulation.xval,FEM_structure.Xmin,FEM_structure.Xmax,formulation.xold1,formulation.xold2, ...
        f0val,df0dx,fval,dfdx,FEM_structure.low,FEM_structure.upp,FEM_structure.a0,FEM_structure.a,FEM_structure.c,FEM_structure.d);
    formulation.xold2 = formulation.xold1;
    formulation.xold1 = formulation.xval;
    formulation.xval  = xmma;
    formulation.x  = xmma;
    change=norm(formulation.xval-formulation.xold1,1);
    [formulation.xPhys,formulation.den_grad] = mna_density(FEM_structure,formulation.x(:)) ;
    %     if FEM_structure.ft == 1
    %         formulation.xPhys=formulation.xval;
    %     elseif FEM_structure.ft == 2
    %         formulation.xPhys = (FEM_structure.H*formulation.xval(:))./FEM_structure.Hs;
    %     end
    %     xPhys(x1>=2500&x4<=4500&z1==min(z1))=0;
    %% penalty correction
    if rem(FEM_structure.outit,20000000)==0
        FEM_structure.Emin=FEM_structure.Emin0/10^(fix(FEM_structure.outit/20));
        FEM_structure.penal=FEM_structure.penal*1.05;
    end
    %% FE-ANALYSIS
    mgcg_selector=false;
    Fac_color=zeros(6*length(formulation.xPhys),1); % TO DO: update?
    for nn=1:length(formulation.xPhys) % TO DO: update?
        Fac_color(6*(nn-1)+(1:6))=(1-formulation.xPhys(nn))*ones(6,1); % TO DO: update?
    end
    Optistruc_affich=0.3;
    Maff=1-Optistruc_affich;
    %color filter
    FACES_shown=FEM_structure.FACES(Fac_color<=Maff,:);
    Fac_color=Fac_color(Fac_color<=Maff);
    FEM_structure.outit=FEM_structure.outit;
    FEM_structure.FACES_shown=FACES_shown;
    [C,dC,perfo,dperfo,GV,dGv,GKSl,dGKSl,FEM_structure]=FEA3D_compliance(FEM_structure,formulation.xPhys,formulation.den_grad,FEM_structure.E0,FEM_structure.penal,formulation.volfrac,mgcg_selector,FEM_structure.nfin,FEM_structure.Engine_master,FEM_structure.VMl,FEM_structure.Aggregation_constant,FEM_structure.toll_minres_U,FEM_structure.toll_minres_TSFC,FEM_structure.toll_minres_GKS);
    
    %     %% FILTERING/MODIFICATION OF SENSITIVITIES §§§ no need !!!? §§§
    %     if FEM_structure.ft == 1
    %         dperfo(:) = FEM_structure.H*(formulation.x(:).*dperfo(:))./FEM_structure.Hs./max(1e-3,formulation.x(:));
    %         dGKSl(:)= FEM_structure.H*(formulation.x(:).*dGKSl(:))./FEM_structure.Hs./max(1e-3,formulation.x(:));
    %         dC(:)= FEM_structure.H*(formulation.x(:).*dC(:))./FEM_structure.Hs./max(1e-3,formulation.x(:));
    %         %     dfandisp(:)=H*(x(:).*dfandisp.')./Hs./max(1e-3,x(:));
    %     elseif FEM_structure.ft == 2
    %         dperfo = FEM_structure.H*(dperfo(:)./FEM_structure.Hs);
    %         dGKSl(:)= FEM_structure.H*(dGKSl(:)./FEM_structure.Hs);
    %         %     dfandisp=H*(dfandisp.'./Hs);
    %         dGv(:) = FEM_structure.H*(dGv(:)./FEM_structure.Hs);
    %         dC(:) = FEM_structure.H*(dC(:)./FEM_structure.Hs);
    %     end
    %% Gather results
    Gc=(C-formulation.admissible_compliance)/formulation.admissible_compliance*100;
    Gp=(perfo-formulation.admissible_DTSFC)/formulation.admissible_DTSFC*100;
    dGp=dperfo/formulation.admissible_DTSFC*100;
    dGc=dC/formulation.admissible_compliance*100;
    Vol=(1+GV)*formulation.volfrac;
    dVol=dGv*formulation.volfrac;
    switch formulation.type
        case 1
            f0val=C;
            fval=100*[GV*100;GKSl;Gp];%;1*(-0.38-FAN_Ax_disp)/0.38
            df0dx=1*reshape(dC,[],1);
            dfdx=100*[dGv(:)'*100;dGKSl(:)';dGp(:)'];%;-1*dfandisp(:)'/0.38
        case 2
            f0val=Vol*100;
            fval=10*[Gc;GKSl;Gp];%;1*(-0.38-FAN_Ax_disp)/0.38
            df0dx=100*reshape(dVol,[],1);
            dfdx=10*[dGc(:)';dGKSl(:)';dGp(:)'];%;-1*dfandisp(:)'/0.38
        case 3
            f0val=perfo*100;
            fval=10*[Gc;GKSl;GV];%;1*(-0.38-FAN_Ax_disp)/0.38
            df0dx=100*reshape(dperfo,[],1);
            dfdx=10*[dGc(:)';dGKSl(:)';dGv(:)'];%;-1*dfandisp(:)'/0.38
    end
    Fac_color=zeros(6*length(formulation.xPhys),1); % TO DO: update?
    for nn=1:length(formulation.xPhys) % TO DO: update?
        Fac_color(6*(nn-1)+(1:6))=(1-formulation.xPhys(nn))*ones(6,1); % TO DO: update?
    end
    Optistruc_affich=0.3;
    Maff=1-Optistruc_affich;
    %color filter
    FACES_shown1=FEM_structure.FACES(Fac_color<=Maff,:);
    Fac_color1=Fac_color(Fac_color<=Maff);
    FACES_shown2=FEM_structure.FACES(Fac_color<=0.1,:);
    Fac_color2=Fac_color(Fac_color<=0.1);
    h = figure(2);clf
    set(h,'Color',[1 1 1]);
    axes1 = axes('Parent',h);
    hold on; h24_patch = patch('Parent',axes1,'Vertices',FEM_structure.COORD,'Faces',FACES_shown1,'FaceVertexCData',Fac_color1,'FaceColor','flat'); axis equal; axis off;
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
    if outit == 498 
        savefig('density_before_transition.fig') ;
    end
    if outit == 598
        savefig('density_after_transition.fig') ;
    end
    %     if print_all
    %         print(['configuration_',num2str(outit,'%03d')],'-djpeg')
    %     end
    %         h = figure(6);clf
    %         set(h,'Color',[1 1 1]);
    %         axes1 = axes('Parent',h);
    %         hold on; h24_patch = patch('Parent',axes1,'Vertices',FEM_structure.COORD,'Faces',FACES_shown2,'FaceVertexCData',Fac_color2,'FaceColor','flat'); axis equal; axis off;
    %         quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
    %         quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
    %         quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
    %         text(1000,0,0,'x')
    %         text(0,1000,0,'y')
    %         text(0,0,1000,'z')
    %         view([27.6 18]);
    %         set(axes1,'CameraViewAngle',3.84174456530832,'DataAspectRatio',[1 1 1],...
    %             'PlotBoxAspectRatio',[434 342.3 684.6]);
    fprintf(' It.:%5i C:%11.4f T.:%11.4f Vol.:%7.3f GKSl.:%7.3f kktn.:%7.3f och.:%7.3f\n',outit,C,perfo ,...
        Vol,GKSl,FEM_structure.kktnorm,error);
    %
    %     figure(2)
    %     hold on
    %     plot(outeriter,perfo,'bo','MarkerFaceColor','b')
    %     plot(outeriter,Vol*100,'ro','MarkerFaceColor','r')
    %     plot(outeriter,(1+GKSl)*VMl,'ko','MarkerFaceColor','k')
    %     title(['Convergence volfrac = ',num2str(Vol*100),', C =',num2str(C),',\sigma_{max} \approx ',num2str((1+GKSl)*VMl),' MPa , iter = ', num2str(outit)])
    %     grid on
    %     xlabel('iter')
    %     legend('\Delta TSFC %','Volume Fraction %','\sigma_{max}')
    %     if print_all
    %         print(['convergence_',num2str(outit,'%03d')],'-djpeg')
    %     end
    %     figure(8)
    %     clf
    %     hold on
    %     Xs=reshape(xPhys(Ic),8*(nfin+1),25*(nfin+1),38*(nfin+1));
    %     p = patch(isosurface(Xg,Yg,Zg,Xs,0.5));
    %     p.FaceColor = 'red';
    %     p.EdgeColor = 'none';
    %     camlight
    %     lighting gouraud
    %     quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k')
    %     quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k')
    %     quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k')
    %     text(1000,0,0,'x')
    %     text(0,1000,0,'y')
    %     text(0,0,1000,'z')
    %     view([27.6 18]);
    %     axis equal
    %     axis off
    outit=FEM_structure.outit;
    save([Path,'convergence_',num2str(FEM_structure.outit,'%03d')],'outit','perfo','GKSl','Vol','C')
    coordo=FEM_structure.COORD;
    save([Path,'configuration_',num2str(FEM_structure.outit,'%03d')],'coordo','FACES_shown')
    clear coordo
    %% The residual vector of the KKT conditions is calculated:
    [residu,FEM_structure.kktnorm,residumax] = ...
        kktcheck(FEM_structure.m,FEM_structure.n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        FEM_structure.xmin,FEM_structure.xmax,df0dx,fval,dfdx,FEM_structure.a0,FEM_structure.a,FEM_structure.c,FEM_structure.d);
    outvector1 = [FEM_structure.outeriter innerit f0val fval(:)'];
    outvector2 = formulation.xval';
    error=abs((fold-f0val))/f0val*100;
    error_old=error;
    %% covergence history
    if outit > 3
        figure(f3)
        subplot(3,2,1)
        scatter(outit,C,4) ; xlabel('iteration') ; ylabel('compliance');
        title('compliance history') ;
        hold on
        subplot(3,2,2)
        scatter(outit,FEM_structure.kktnorm,4) ; xlabel('iteration') ; ylabel('KKT norm');
        title('KKT norm history') ;
        hold on
        subplot(3,2,3)
        scatter(outit,Vol,4) ; xlabel('iteration') ; ylabel('mass constraint');
        title('mass constraint history') ;
        hold on
        subplot(3,2,4)
        scatter(outit,perfo,4) ; xlabel('iteration') ; ylabel('TSFC constraint');
        title('TSFC history') ;
        hold on
        subplot(3,2,5)
        scatter(outit,GKSl,4) ; xlabel('iteration') ; ylabel('stress aggregation');
        title('stress constraint history') ;
        hold on
    end
    if outit == 498
        savefig('history_before_transition.fig') ;
    end
    if outit == 598
        savefig('history_after_transition.fig') ;
    end
end
kktnorm=FEM_structure.kktnorm;
x=formulation.xval; % §§§§§
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
