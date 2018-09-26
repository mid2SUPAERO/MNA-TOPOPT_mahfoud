function [perfo,FAN_Ax_disp,TLs_intersection]=perfo_func(x,T)
tic
[K]=Stiffness_mat_assembly(x,T);
F=F(DOFs_MAP(:,1)~=P2(1)&DOFs_MAP(:,1)~=P2d(1),:);
Final_DOF_MAP=DOFs_MAP(DOFs_MAP(:,1)~=P2(1)&DOFs_MAP(:,1)~=P2d(1),:);
%Displacement evaluation
U=K\F;
%Disp-tip clearance
U_out=U(1:Number_of_retained_DOF,:);
Utip=Recovery_Matrix*U_out;
load('gamma_structure')
load('U_tilder')
perfo=0;
for s=1:20
tip_vector=Gamma(s).mat*(U_static+Utip);
R=rms(tip_vector);
perfo=perfo+Gamma(s).lambda*R;
end
perfo=perfo(1);
FAN_Ax_disp=U_FAN_static+FAN_Center_Recovery_vector*U_out;
FAN_Ax_disp=FAN_Ax_disp(1);
%
TLs_intersection=[P1(3),P1(4),1;P2(3),P2(4),1;P1d(3),P1d(4),1]\[-P1(2);-P2(2);-P1d(2)];
TLs_intersection=-TLs_intersection(3);
toc