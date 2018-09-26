function [perfo,FAN_Ax_disp,TLs_intersection,grad_perfo,grad_FAN_Ax_disp]=perfo_func_sensitivity(x,T)
tic
load('substructure_results.mat')
[K,F,DOFs_MAP,P1,P1d,P2,P2d,Number_of_retained_DOF]=Stiffness_mat_assembly(x,T);
eps_x=0.1;
eps_T=0.1;
%stiffness_sensitiity
[Kxp]=Stiffness_mat_assembly(x+eps_x,T);
[Kxm]=Stiffness_mat_assembly(x-eps_x,T);
[KTm]=Stiffness_mat_assembly(x,T-eps_T);
[KTp]=Stiffness_mat_assembly(x,T+eps_T);
%
DK_DX=(Kxp-Kxm)/2/eps_x;
DK_DT=(KTp-KTm)/2/eps_T;
F=F(DOFs_MAP(:,1)~=P2(1)&DOFs_MAP(:,1)~=P2d(1),:);
Final_DOF_MAP=DOFs_MAP(DOFs_MAP(:,1)~=P2(1)&DOFs_MAP(:,1)~=P2d(1),:);
%Displacement evaluation
F=[F,[Recovery_Matrix.';zeros(size(F,1)-Number_of_retained_DOF,size(Recovery_Matrix,1))],[FAN_Center_Recovery_vector.';zeros(size(F,1)-Number_of_retained_DOF,size(FAN_Center_Recovery_vector,1))]];
Res=K\F;
U=Res(:,1:3);
H=Res(:,4:end-1);
H_f=Res(:,end);
%Disp-tip clearance
U_out=U(1:Number_of_retained_DOF,:);
Utip=Recovery_Matrix*U_out;
load('gamma_structure')
load('U_tilder')
perfo=0;
psi_tild=0;

for s=1:20
tip_vector=Gamma(s).mat*(U_static+Utip);
R=rms(tip_vector);
perfo=perfo+Gamma(s).lambda*R;
N=length(tip_vector);
R=ones(N,1)*R;
psi_tild=psi_tild-Gamma(s).lambda/N*H*Gamma(s).mat.'*(tip_vector./R);
end
perfo=perfo(1);
dperfo_dx=-psi_tild.'*DK_DX*U;
dperfo_dx=dperfo_dx(1,1);
dperfo_dT=-psi_tild.'*DK_DT*U;
dperfo_dT=dperfo_dT(1,1);
grad_perfo=[dperfo_dx;dperfo_dT];
FAN_Ax_disp=U_FAN_static+FAN_Center_Recovery_vector*U_out;
FAN_Ax_disp=FAN_Ax_disp(1);
dFAN_Ax_disp_dx=-H_f.'*DK_DX*U;
dFAN_Ax_disp_dx=dFAN_Ax_disp_dx(1);
dFAN_Ax_disp_dT=-H_f.'*DK_DT*U;

dFAN_Ax_disp_dT=dFAN_Ax_disp_dT(1);
grad_FAN_Ax_disp=[dFAN_Ax_disp_dx;dFAN_Ax_disp_dT];
%
TLs_intersection=[P1(3),P1(4),1;P2(3),P2(4),1;P1d(3),P1d(4),1]\[-P1(2);-P2(2);-P1d(2)];
TLs_intersection=-TLs_intersection(3);
toc