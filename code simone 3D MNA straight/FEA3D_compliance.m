function [c,dc,perfo,dperfo,GV,dGv,GKSl,dGKSl,FEM_structure]=FEA3D_compliance(FEM_structure,xPhys,den_grad,E0,penal,volfrac,mgcg_selector,nl,Engine_master,VMl,Penalty,toll_minres_U,toll_minres_TSFC,toll_minres_GKs) % TO DO: update?
% inputs
free_dofs=FEM_structure.free_dofs;
K_def=FEM_structure.K_def;
A2=FEM_structure.A2;
U=FEM_structure.U;
Sol=FEM_structure.Sol;
psi_tild_sol=FEM_structure.psi_tild_sol;
psi_tild=FEM_structure.psi_tild;
if mgcg_selector
    Null=FEM_structure.Null;
end
%% FE-ANALYSIS
E=FEM_structure.Emin+xPhys(:).^penal*(E0-FEM_structure.Emin);
E=E(FEM_structure.el_ij);
sK=FEM_structure.k_ij.*E;
% X=xPhys(FEM_structure.elm_ij);
% Mk=FEM_structure.M_ij.*X;

if mgcg_selector
    K=cell(nl,1);
    K{1,1} = sparse(FEM_structure.I,FEM_structure.J,sK,size(K_def,1),size(K_def,2)); K{1,1}=K{1,1}+triu(K{1,1}.',1);
    K{1,1}=K_def+K{1,1};
    K{1,1} = Null'*K{1,1}*Null - (Null-speye(size(K{1,1})));
else
    if ~Engine_master
        K = sparse(FEM_structure.I,FEM_structure.J,sK,size(K_def,1),size(K_def,2)); K=K+triu(K.',1); K=(K+K.')/2;
        K=K_def+K;
    else
        K = sparse(FEM_structure.I,FEM_structure.J,sK,size(K_def,1),size(K_def,2)); K=K+triu(K.',1); K=(K+K.')/2;
        K=FEM_structure.K_Engine+FEM_structure.Pr_EM'*K*FEM_structure.Pr_EM;
    end
end
if ~Engine_master
    F=(FEM_structure.LOAD_MATRIX_def(:,1));
else
    F=FEM_structure.LOAD_MATRIX_EM(:,1);
end

if mgcg_selector
    % Kfac=chol(K(free_dofs,free_dofs));
    % Res=Kfac\(Kfac.'\F(free_dofs,:));
    for l = 1:nl-1
        %     Put=[Pu{l,1},zeros(size(Pu{l,1},1),sum(DOFs_MAP(:,1)==0));zeros(sum(DOFs_MAP(:,1)==0),size(Pu{l,1},2)),speye(sum(DOFs_MAP(:,1)==0))];
        %     Put{l,1}=Pu{l,1}(1:(end-sum(DOFs_MAP(:,1)==0)),1:(end-sum(DOFs_MAP(:,1)==0)));
        K{l+1,1} = FEM_structure.Pu{l,1}'*(K{l,1}*FEM_structure.Pu{l,1});
        %     K_tild{l+1,1} = Pu{l,1}'*(K_tild{l,1}*Pu{l,1});
    end
    Lfac = chol((K{nl,1}),'lower'); Ufac = Lfac';
    cgtol=1e-10;cgmax=100;
    tic
    [~,~,U] = mgcg(K,F,U,Lfac,Ufac,FEM_structure.Pu,nl,1,cgtol,cgmax);
else
    if length(free_dofs)<1e4||toll_minres_U==0
        if ~Engine_master
            Res=K(free_dofs,free_dofs)\F(free_dofs,:);U(free_dofs)=Res(:,1);
            Fdef=K_def*U;
            c=U'*K*U-U'*Fdef;
        else
            Res=K(FEM_structure.free_dofs_engine_master,FEM_structure.free_dofs_engine_master)\F(FEM_structure.free_dofs_engine_master,:);Sol(FEM_structure.free_dofs_engine_master)=Res;
            U=FEM_structure.Pr_EM*Sol;
            Fdef=FEM_structure.K_Engine*Sol;
            c=Sol'*K*Sol-Sol'*Fdef;
        end
    else
        if ~Engine_master
            alpha = (max(sum(abs(K(free_dofs,free_dofs)),2)./diag(K(free_dofs,free_dofs)))-2)/100;
%                     alpha=0.1;
            L1 = sparse(ichol(K(free_dofs,free_dofs), struct('type','ict','droptol',1e-3,'diagcomp',alpha/10)));
            maxit=1000;
            [Res,fl1,rr1,it1,rv1] = minres(K(free_dofs,free_dofs),F(free_dofs,:),1e-10,maxit,L1,L1',U(free_dofs));
            if fl1~=0
                disp(['pcg didn''t converge to the desired tolerance relres =',num2str(rr1)])
            end
            %         L = ichol(K(free_dofs,free_dofs),struct('michol','on'));
            %         [Res,fl2,rr2,it2,rv2] = cgs(K(free_dofs,free_dofs),F(free_dofs,:),1e-8,10000);
            U(free_dofs)=Res(:,1);
            Fdef=K_def*U;
            c=U'*K*U-U'*Fdef;
        else
            alpha = (max(sum(abs(K(FEM_structure.free_dofs_engine_master,FEM_structure.free_dofs_engine_master)),2)./diag(K(FEM_structure.free_dofs_engine_master,FEM_structure.free_dofs_engine_master)))-2)/100;
%                     alpha=0.1;
            tic
            r = symrcm(K(FEM_structure.free_dofs_engine_master,FEM_structure.free_dofs_engine_master));
            L1 =sparse( ichol(K(FEM_structure.free_dofs_engine_master(r),FEM_structure.free_dofs_engine_master(r)), struct('type','ict','droptol',1e-3,'diagcomp',alpha/10)));
%             L1=diag(diag(K(FEM_structure.free_dofs_engine_master,FEM_structure.free_dofs_engine_master)));
            maxit=1000;
            [Res,fl1,rr1,it1,rv1] = minres(K(FEM_structure.free_dofs_engine_master(r),FEM_structure.free_dofs_engine_master(r)),F(FEM_structure.free_dofs_engine_master(r)),toll_minres_U,maxit,L1,L1',Sol(FEM_structure.free_dofs_engine_master(r)));
            toc
            if fl1~=0
                disp(['pcg didn''t converge to the desired tolerance relres =',num2str(rr1)])
            end
            %         L = ichol(K(free_dofs,free_dofs),struct('michol','on'));
            %         [Res,fl2,rr2,it2,rv2] = cgs(K(free_dofs,free_dofs),F(free_dofs,:),1e-8,10000);
            Sol(FEM_structure.free_dofs_engine_master(r))=Res;
            U=FEM_structure.Pr_EM*Sol;
            Fdef=FEM_structure.K_Engine*Sol;
            c=Sol'*K*Sol-Sol'*Fdef;
        end
%         %% Test pcg mires bcg
%                 figure(10)
%                 semilogy(0:max(it1,length(rv1)-1),rv1./norm(F(FEM_structure.free_dofs_engine_master,:)),'r.');
%         %         semilogy(0:10000,rv2/norm(F(free_dofs,:)),'-o');
%                 xlabel('Iteration number');
%                 ylabel('Relative residual');
    end
end

% L = ichol(K(free_dofs,free_dofs),struct('michol','on'));
% [Res(:,1),fl2,rr2,it2,rv2] = cgs(K(free_dofs,free_dofs),F(free_dofs,:),1e-8,100);
% toc
xPhys=xPhys(:);
CdK=(penal*(E0-FEM_structure.Emin)*(xPhys(:)).^(penal-1)); % TO DO: update?
dK_dxU=reshape(FEM_structure.K0c*U,24,[])*spdiags(CdK,0,length(CdK),length(CdK));
dK_dxU=dK_dxU(:);
dK_dxU=sparse(reshape(FEM_structure.DESIGN_ZONE_ELEMENT_DOFS',[],1)...
    ,reshape(repmat(1:size(FEM_structure.DESIGN_ZONE_ELEMENT_DOFS,1),24,1),[],1),dK_dxU);
dK_dxU = dK_dxU*den_grad ;
%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
if ~Engine_master
    Utip=FEM_structure.Recovery_Matrix_def*U;
    perfo=0;
    DTSFC_DU=0;
    for s=1:20
        tip_vector=sparse(FEM_structure.Gamma(s).mat)*(FEM_structure.U_static(:,1)+Utip);
        R=rms(tip_vector);
        perfo=(perfo+FEM_structure.Gamma(s).lambda*R);
        N=length(tip_vector);
        R=ones(N,1)*R;
        DTSFC_DU=DTSFC_DU+FEM_structure.RG{s}*(tip_vector./R);
    end
else
    Utip=FEM_structure.Recovery_Matrix_EM*Sol;
    perfo=0;
    DTSFC_DU=0;
    for s=1:20
        tip_vector=sparse(FEM_structure.Gamma(s).mat)*(FEM_structure.U_static(:,1)+Utip);
        R=rms(tip_vector);
        perfo=(perfo+FEM_structure.Gamma(s).lambda*R);
        N=length(tip_vector);
        R=ones(N,1)*R;
        DTSFC_DU=DTSFC_DU+FEM_structure.RG_EM{s}*(tip_vector./R);
    end
end
% RES2=K(free_dofs,free_dofs).'\[DTSFC_DU(free_dofs),FAN_Center_Recovery_vector_def(free_dofs).'];
micro_stress=FEM_structure.DB*U;
micro_stress=reshape(micro_stress,6,[]);
compl_stress=[micro_stress(1,:)-micro_stress(2,:);micro_stress(2,:)-micro_stress(3,:);micro_stress(3,:)-micro_stress(1,:);3*micro_stress(4,:);3*micro_stress(5,:);3*micro_stress(6,:)];
VM2=sum(compl_stress.*micro_stress);
VM=sqrt(VM2);
RO=repmat(xPhys(:)',8,1);RO=reshape(RO,size(VM2,1),size(VM2,2));
% % VMl=25;
gj=RO.*(VM/VMl-1);
% % Penalty=10;
Ngp=length(RO);
[gmax]=max(gj);
GKSl=(gmax-log(Ngp)/Penalty+1/Penalty*log(sum(exp(Penalty*(gj-gmax)))));
mccoef=([1 -0.5 -0.5 0 0 0;-0.5 1 -0.5 0 0 0;-0.5 -0.5 1 0 0 0;0 0 0 3 0 0;0 0 0 0 3 0;0 0 0 0 0 3]*micro_stress).*repmat(RO./VM.*exp(Penalty*(gj-gmax)),6,1)/VMl/(sum(exp(Penalty*(gj-gmax))));
if ~Engine_master
    Tau=mccoef(:).'*FEM_structure.DB;
else
    Tau=mccoef(:).'*FEM_structure.DB*FEM_structure.Pr_EM;
end
% RES2=Kfac\(Kfac.'\[DTSFC_DU(free_dofs),Tau(free_dofs)']);

% RES2=K(free_dofs,free_dofs)\[DTSFC_DU(free_dofs),Tau(free_dofs)'];
if ~mgcg_selector
    if length(free_dofs)<1e4||toll_minres_U==0
        if ~Engine_master
            RES2=K(free_dofs,free_dofs)\[DTSFC_DU(free_dofs),Tau(free_dofs)',2*Fdef(free_dofs)];psi_tild(free_dofs,:)=RES2;
%             RES2=K(free_dofs,free_dofs)\[DTSFC_DU(free_dofs)];
        else
            RES2=K(FEM_structure.free_dofs_engine_master,FEM_structure.free_dofs_engine_master)\[DTSFC_DU(FEM_structure.free_dofs_engine_master,:),Tau(FEM_structure.free_dofs_engine_master)',2*Fdef(FEM_structure.free_dofs_engine_master)];psi_tild_sol(FEM_structure.free_dofs_engine_master,:)=RES2;
            psi_tild=FEM_structure.Pr_EM*psi_tild_sol;
        end
    else
        %          alpha=0.1;
        %         L1 = ichol(K(free_dofs,free_dofs), struct('type','ict','droptol',1e-3,'diagcomp',alpha));
        if ~Engine_master
            maxit=1000;
            [RES2(:,1),fl1,rr1,it1,rv1] = minres(K(free_dofs,free_dofs),DTSFC_DU(free_dofs),1e-10,maxit,L1,L1',psi_tild(free_dofs,1));psi_tild(free_dofs,1)=RES2(:,1);
            if fl1~=0
                disp(['pcg didn''t converge to the desired tolerance relres =',num2str(rr1)])
            end
            [RES2(:,2),fl1,rr1,it1,rv1] = minres(K(free_dofs,free_dofs),Tau(free_dofs)',1e-10,maxit,L1,L1',psi_tild(free_dofs,2));psi_tild(free_dofs,2)=RES2(:,2);
            if fl1~=0
                disp(['pcg didn''t converge to the desired tolerance relres =',num2str(rr1)])
            end
            [RES2(:,3),fl1,rr1,it1,rv1] = minres(K(free_dofs,free_dofs),2*Fdef(free_dofs),1e-10,maxit,L1,L1',psi_tild(free_dofs,3));psi_tild(free_dofs,3)=RES2(:,3);
            if fl1~=0
                disp(['pcg didn''t converge to the desired tolerance relres =',num2str(rr1)])
            end
        else
            maxit=1000;
%             DTSFC_DU=FEM_structure.Pr_EM'*DTSFC_DU;
            [RES2(:,1),fl1,rr1,it1,rv1] = minres(K(FEM_structure.free_dofs_engine_master(r),FEM_structure.free_dofs_engine_master(r)),DTSFC_DU(FEM_structure.free_dofs_engine_master(r)),toll_minres_TSFC,maxit,L1,L1',psi_tild_sol(FEM_structure.free_dofs_engine_master(r),1));psi_tild_sol(FEM_structure.free_dofs_engine_master(r),1)=RES2(:,1);
            if fl1~=0
                disp(['pcg didn''t converge to the desired tolerance relres =',num2str(rr1)])
            end
            
            [RES2(:,2),fl1,rr1,it1,rv1] = minres(K(FEM_structure.free_dofs_engine_master(r),FEM_structure.free_dofs_engine_master(r)),Tau(FEM_structure.free_dofs_engine_master(r))' ,toll_minres_GKs,maxit,L1,L1',psi_tild_sol(FEM_structure.free_dofs_engine_master(r),2));psi_tild_sol(FEM_structure.free_dofs_engine_master(r),2)=RES2(:,2);
            if fl1~=0
                disp(['pcg didn''t converge to the desired tolerance relres =',num2str(rr1)])
            end
            [RES2(:,3),fl1,rr1,it1,rv1] = minres(K(FEM_structure.free_dofs_engine_master(r),FEM_structure.free_dofs_engine_master(r)),2*Fdef(FEM_structure.free_dofs_engine_master(r)),toll_minres_GKs,maxit,L1,L1',psi_tild_sol(FEM_structure.free_dofs_engine_master(r),3));psi_tild_sol(FEM_structure.free_dofs_engine_master(r),3)=RES2(:,3);
            if fl1~=0
                disp(['pcg didn''t converge to the desired tolerance relres =',num2str(rr1)])
            end
            psi_tild=FEM_structure.Pr_EM*psi_tild_sol;
        end
        %         [RES2,fl2,rr2,it2,rv2] = cgs(K(free_dofs,free_dofs),DTSFC_DU(free_dofs),1e-8,100);psi_tild(free_dofs,1)=RES2;
    end
else
    
    [~,~,psi_tild(:,1)] = minres(K,DTSFC_DU,FEM_structure.psi_tild(:,1),Lfac,Ufac,FEM_structure.Pu,nl,1,cgtol,cgmax);
end
% [~,~,psi_tild(:,2)] = mgcg(K,[Tau'],FEM_structure.psi_tild(:,2),Lfac,Ufac,FEM_structure.Pu,nl,1,cgtol,cgmax);



% FAN_Ax_disp=U_FAN_static(1)+FAN_Center_Recovery_vector_def*U;
dperfo=(-psi_tild(:,1).'*dK_dxU); % TO DO: update?
dGKSl_dro=+1/(sum(exp(Penalty*(gj-gmax))))*sum(reshape((exp(Penalty*(gj-gmax)).*(VM/VMl-1))',8,[])); 
dGKSl_dro = dGKSl_dro*den_grad ; % dGKSl instead; derivatives are wrt design variables that are no more densities
dGKSl=-psi_tild(:,2).'*dK_dxU+dGKSl_dro;
dc=(psi_tild(:,3)-U)'*dK_dxU;
% dfandisp=-psi_tild(:,2).'*dK_dxU;
dv = A2/2;
dGv=(dv(:)'/(volfrac*sum(A2/2)))*den_grad; % TO DO: update?
GV=(sum(xPhys(:).*A2/2) - volfrac*sum(A2/2))/(volfrac*sum(A2/2));
if FEM_structure.stress_plot==1
    VMp=(VM(:)/VMl-1).*RO(:);
    VMn=(FEM_structure.Ngpn(1:6:end,1:6:end)'*FEM_structure.Ngpn(1:6:end,1:6:end))\(FEM_structure.Ngpn(1:6:end,1:6:end)'*VMp(:));
    ColorOdsR = VMn;
     if ~FEM_structure.save_res
    h = figure(6);clf
    set(h,'Color',[1 1 1]);
    axes1 = axes('Parent',h);
    hold on;  patch('Parent',axes1,'Vertices',FEM_structure.COORD,'Faces',FEM_structure.FACES_shown,'FaceVertexCData',ColorOdsR,'FaceColor','interp','EdgeColor','none'); axis equal; axis off; drawnow;
    colormap jet;
    axis equal
    axis off
    s1=quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k');
    s1=quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k');
    s1=quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k');
    text(1000,0,0,'x')
    text(0,1000,0,'y')
    text(0,0,1000,'z')
    view([27.6 18]);
    % Set the remaining axes properties
    set(axes1,'CameraViewAngle',3.84174456530832,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatio',[434 342.3 684.6]);
    C=colorbar;
    caxis([-1 0.5])
    C.Limits=[-1 0.5];
    C.Label.String = 'x(\sigma_{VM}/\sigma_{lim}-1) ';
    C.FontWeight='bold';
    C.Location='southoutside';
    drawnow;
    title(['Von Mises stress at iteration ',num2str(FEM_structure.outit)])
    if FEM_structure.print_all
        print(['VMstress_',num2str(FEM_structure.outit,'%03d')],'-dpng')
    end
    else
        save(['VMstress_',num2str(FEM_structure.outit,'%03d')],'ColorOdsR')
    end
    
end
FEM_structure.U=U;
FEM_structure.psi_tild=psi_tild;
FEM_structure.Sol=Sol;
FEM_structure.psi_tild_sol=psi_tild_sol;
