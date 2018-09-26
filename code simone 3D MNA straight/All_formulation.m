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
%
print_all=~true;
filter=true;
save_res= ~true;
FEM_structure.save_res=save_res;
FEM_structure.stress_plot=0;
FEM_structure.print_all=print_all;
% Nthr=maxNumCompThreads(8);
FEM_structure.Engine_master=true;
plot_problem=true;
%the point at the same x have the same 1 and 3 diplacement + linear
%% check if there are Lagrange multiplicator and eliminate orphan nodes
[FEM_structure] = No_lag(FEM_structure);
%% refinement
FEM_structure.nfin=0;
tic
[FEM_structure]=mesh_refinement_3D(FEM_structure,FEM_structure.nfin);
toc
%% Engine dofs eliminaion through fem shape function interpolation
%the point at the same x have the same 1 and 3 diplacement + linear
tic
[FEM_structure]= engine_DOFs_eliminationpar(FEM_structure);
toc
%% plot the engine and DZ mesh together
El_numb=size(FEM_structure.ELEMENT,1);
FEM_structure.FACES=zeros(6*El_numb,4);
FID=[1 2 3 4;5 6 7 8;2 3 7 6;3 4 8 7;4 1 5 8;1 2 6 5];
for nn=1:El_numb
    FEM_structure.FACES(6*(nn-1)+(1:6),1)=FEM_structure.ELEMENT(nn,FID(:,1));
    FEM_structure.FACES(6*(nn-1)+(1:6),2)=FEM_structure.ELEMENT(nn,FID(:,2));
    FEM_structure.FACES(6*(nn-1)+(1:6),3)=FEM_structure.ELEMENT(nn,FID(:,3));
    FEM_structure.FACES(6*(nn-1)+(1:6),4)=FEM_structure.ELEMENT(nn,FID(:,4));
end
% FEM_structure.FACES=FACES;
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
    hold on; h24_patch = patch('Parent',axes1,'Vertices',FEM_structure.COORD,'Faces',FEM_structure.FACES,'FaceVertexCData',[.5 .1 .1],'FaceColor','flat'); axis equal; axis off;
    hold on; h25_patch = patch('Parent',axes1,'Vertices',FEM_structure.intercoord(:,2:4),'Faces',FACE_intE(:,2:5),'FaceVertexCData',[.6 .5 .3],'FaceColor','flat'); axis equal; axis off;
    s1=quiver3(0,0,0,1000,0,0,'LineWidth',3,'Color','k');
    s2=quiver3(0,0,0,0,1000,0,'LineWidth',3,'Color','k');
    s3=quiver3(0,0,0,0,0,1000,'LineWidth',3,'Color','k');
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
FEM_structure.E0=210000;
[FEM_structure]= Stiffness_components3D_compliance(FEM_structure,FEM_structure.E0);
toc
%% INITIALIZE ITERATION

FEM_structure.penal=3; FEM_structure.ft=2;



% xPhys(x1>=3000&x4<=5000&z1<=min(z1)+dz)=0;
FEM_structure.loop = 0;
FEM_structure.change = 1;
FEM_structure.m = 3;
FEM_structure.n = size(FEM_structure.DESIGN_ZONE_ELEMENT_NODES,1);
FEM_structure.epsimin = 0.0000001;
FEM_structure.nX = 5 ; % number of components in X
FEM_structure.nY = 5 ; % number of components in Y
FEM_structure.nZ = 5 ; % number of components in Z
FEM_structure.n_var_c = 9 ; % number of variables per components
FEM_structure.n_var = FEM_structure.n_var_c*FEM_structure.nZ*...
    FEM_structure.nY*FEM_structure.nX ; % total number of components
FEM_structure.eeen    = ones(FEM_structure.n_var,1);
FEM_structure.eeem    = ones(FEM_structure.m,1);
FEM_structure.zeron   = zeros(FEM_structure.n_var,1);
FEM_structure.zerom   = zeros(FEM_structure.m,1);
FEM_structure.limits = zeros(1,6) ; % limites of design domain
FEM_structure.limits(1) = min(FEM_structure.COORD(:,1)) ;
FEM_structure.limits(2) = max(FEM_structure.COORD(:,1)) ;
FEM_structure.limits(3) = min(FEM_structure.COORD(:,2)) ;
FEM_structure.limits(4) = max(FEM_structure.COORD(:,2)) ;
FEM_structure.limits(5) = min(FEM_structure.COORD(:,3)) ;
FEM_structure.limits(6) = max(FEM_structure.COORD(:,3)) ;
FEM_structure.Xmin = [FEM_structure.limits(1);FEM_structure.limits(3);...
    FEM_structure.limits(5);0;0;0;-pi;-pi;-pi] ;
FEM_structure.Xmin = repmat(FEM_structure.Xmin,FEM_structure.n_var/...
    FEM_structure.n_var_c,1) ;
FEM_structure.Xmax = [FEM_structure.limits(2);FEM_structure.limits(4);...
    FEM_structure.limits(6);2*(FEM_structure.limits(2)-...
    FEM_structure.limits(1));2*(FEM_structure.limits(4)-...
    FEM_structure.limits(3));2*(FEM_structure.limits(6)-...
    FEM_structure.limits(5));pi;pi;pi] ;
FEM_structure.Xmax = repmat(FEM_structure.Xmax,FEM_structure.n_var/...
    FEM_structure.n_var_c,1) ;
FEM_structure.xmin    = FEM_structure.zeron;
FEM_structure.xmax    = FEM_structure.eeen;
FEM_structure.low     = FEM_structure.Xmin;
FEM_structure.upp     = FEM_structure.Xmax;
FEM_structure.c       = 1000*FEM_structure.eeem;
FEM_structure.d       = FEM_structure.eeem;
FEM_structure.a0      = 1;
FEM_structure.a       = FEM_structure.zerom;
FEM_structure.outeriter = 0;
FEM_structure.maxoutit  = 600;
FEM_structure.kkttol  = 1e-2;
%%%% The iterations start:
FEM_structure.kktnorm = FEM_structure.kkttol+1000;
FEM_structure.outit = 0;
%% Preparing Filter
% FEM_structure.rmin=1;
% tic
% [FEM_structure.H,FEM_structure.Hs]=preparing_Filter3Dpar(FEM_structure,FEM_structure.rmin);
% toc
%% Optimization
%start by only stress-constrained problem
formulation(1).type=1;
formulation(1).volfrac=0.1;
formulation(1).admissible_DTSFC=0.15;
formulation(1).admissible_compliance=1e7;
% formulation(1).x = ones(size(FEM_structure.DESIGN_ZONE_ELEMENT_NODES,1),1);
% x = volfrac*ones(size(FEM_structure.DESIGN_ZONE_ELEMENT_NODES,1),1);
formulation(1).x = initial_design(FEM_structure) ; % MNA variables
x = initial_design(FEM_structure) ;
[formulation(1).xPhys,formulation(1).den_grad] = mna_density(FEM_structure,formulation(1).x(:)) ;  % xPhys = density field
formulation(1).xval    = formulation(1).x(:); 
formulation(1).xold1   = formulation(1).xval;
formulation(1).xold2   = formulation(1).xval;
% for i=1:3
%     formulation(i).type=i;
%     formulation(i).volfrac=1;
%     formulation(i).admissible_DTSFC=10;
%     formulation(i).admissible_compliance=1e7;
%     formulation(i).x = ones(size(FEM_structure.DESIGN_ZONE_ELEMENT_NODES,1),1);
%     % x = volfrac*ones(size(FEM_structure.DESIGN_ZONE_ELEMENT_NODES,1),1);
%     formulation(i).xPhys = formulation(i).x(:);
%     formulation(i).xval    = formulation(i).xPhys(:);
%     formulation(i).xold1   = formulation(i).xval;
%     formulation(i).xold2   = formulation(i).xval;
% end
i=1;
[C(i),perfo(i),Vol(i),x(:,i),kktnorm(i)]=MMA_optimization(FEM_structure,formulation(i));

% % figure
% % scatter3(C,perfo,Vol,'k','fill')
% % xlabel('C')
% % ylabel('T')
% % zlabel('V')
% FEM_structure.maxoutit  = 200;
% %continue with intermediate solution (2 per side)
% nl=3;%number of level
% %loop over the goal
% formulation=struct([]);
% for i=1:3
%     constraint_set=setdiff(1:3,i);
%     %loop over the first constraint
%     for j=1:2
%         active_constraint=constraint_set(j);
%         switch active_constraint
%             case 1
%                 C0=linspace(C(i),C(active_constraint),nl+2);
%                 C0=setdiff(C0,[C(i),C(active_constraint)]);
%             case 2
%                 V0=linspace(Vol(i),Vol(active_constraint),nl+2);
%                 V0=setdiff(V0,[Vol(i),Vol(active_constraint)]);
%             case 3
%                 T0=linspace(perfo(i),perfo(active_constraint),nl+2);
%                 T0=setdiff(T0,[perfo(i),perfo(active_constraint)]);
%         end
%         %loop over the pareto discretization
%         for ndisc=1:nl
%             formulation(i,j,ndisc).type=i;
%             switch active_constraint
%                 case 1
%                     formulation(i,j,ndisc).volfrac=1;
%                     formulation(i,j,ndisc).admissible_DTSFC=10;
%                     formulation(i,j,ndisc).admissible_compliance=C0(nl);
%                 case 2
%                     formulation(i,j,ndisc).volfrac=V0(nl);
%                     formulation(i,j,ndisc).admissible_DTSFC=10;
%                     formulation(i,j,ndisc).admissible_compliance=1e7;
%                 case 3
%                     formulation(i,j,ndisc).volfrac=1;
%                     formulation(i,j,ndisc).admissible_DTSFC=T0(nl);
%                     formulation(i,j,ndisc).admissible_compliance=1e7;
%             end
%             formulation(i,j,ndisc).xPhys = x(:,active_constraint);
%             formulation(i,j,ndisc).xval    = x(:,active_constraint);
%             formulation(i,j,ndisc).xold1   = x(:,active_constraint);
%             formulation(i,j,ndisc).xold2   = x(:,active_constraint);
%
%         end
%     end
% end
% formulation=formulation(:);
% parfor i=4:(length(formulation)+3)
%    [C(i),perfo(i),Vol(i),x(:,i),kktnorm(i)]=MMA_optimization(FEM_structure,formulation(i-3));
% end
save('final_result','C','perfo','Vol','x','kktnorm','formulation')
% figure
% scatter3(C,perfo,Vol,'k','fill')
% xlabel('C')
% ylabel('T')
% zlabel('V')
% hold on
% quiver3(C(2),perfo(2),Vol(2),C(4)-C(2),perfo(4)-perfo(2),Vol(4)-Vol(2),0)
% quiver3(C(1),perfo(1),Vol(1),C(5)-C(1),perfo(5)-perfo(1),Vol(5)-Vol(1),0)
% quiver3(C(1),perfo(1),Vol(1),C(6)-C(1),perfo(6)-perfo(1),Vol(6)-Vol(1),0)
% quiver3(C(3),perfo(3),Vol(3),C(7)-C(3),perfo(7)-perfo(3),Vol(7)-Vol(3),0)
% quiver3(C(3),perfo(3),Vol(3),C(8)-C(3),perfo(8)-perfo(3),Vol(8)-Vol(3),0)
% quiver3(C(2),perfo(2),Vol(2),C(9)-C(2),perfo(9)-perfo(2),Vol(9)-Vol(2),0)
% axis([min(C) max(C) min(perfo) max(perfo) min(Vol) max(Vol)])