function [K,F,DOFs_MAP,P1,P1d,P2,P2d,Number_of_retained_DOF]=Stiffness_mat_assembly(x,T)
disp('Stiffness Matrix Assembly')
tic
load('substructure_results.mat')
load('lin_pert_sub.dat.mat')
load('connectivity_mat.mat')
Radius=norm([668.671021,121.383797]);
P1=[(max(intercoord(:,1)+1)) 2.451298096E+03,-Radius*sind(T),Radius*cosd(T)];
P2=[(max(intercoord(:,1)+2)) x,-66.380302,910.036987];
P1d=[(max(intercoord(:,1)+3)) 2.451298096E+03,Radius*sind(T),Radius*cosd(T)];
P2d=[(max(intercoord(:,1)+4)) x,66.380302,910.036987];
P=[P1;P2;P1d;P2d];
COORD=[intercoord;P];
%Calculation of Beam and Rod length
Beam_Length=sqrt((COORD(LOCAL_BEAM(:,3),2)-COORD(LOCAL_BEAM(:,2),2)).^2+(COORD(LOCAL_BEAM(:,3),3)-COORD(LOCAL_BEAM(:,2),3)).^2+(COORD(LOCAL_BEAM(:,3),4)-COORD(LOCAL_BEAM(:,2),4)).^2);
Rod_Length=sqrt((COORD(LOCAL_ROD(:,3),2)-COORD(LOCAL_ROD(:,2),2)).^2+(COORD(LOCAL_ROD(:,3),3)-COORD(LOCAL_ROD(:,2),3)).^2+(COORD(LOCAL_ROD(:,3),4)-COORD(LOCAL_ROD(:,2),4)).^2);
phase_square_difference=(atan2(COORD(LOCAL_BEAM(:,2),4),COORD(LOCAL_BEAM(:,2),3))-atan2(COORD(LOCAL_BEAM(:,3),4),COORD(LOCAL_BEAM(:,3),3))).^2;
%Approximation of Heaviside function for Young module evaluation
betha=50;threshold=0.015;
H=@(x) 1-1./(1+exp(-2*betha*(x-threshold)));
%Young Module Evaluation
Emax=1e9;
E=Emax*H(phase_square_difference);
%Rod Young Module
Erod=210000;
%Beam Elementary matrix and Local to Global transformation matrix
%Euler angle evaluation
theta=atan2(COORD(LOCAL_BEAM(:,3),3)-COORD(LOCAL_BEAM(:,2),3),COORD(LOCAL_BEAM(:,3),2)-COORD(LOCAL_BEAM(:,2),2));
psi=atan2(COORD(LOCAL_BEAM(:,3),4)-COORD(LOCAL_BEAM(:,2),4),sqrt((COORD(LOCAL_BEAM(:,3),2)-COORD(LOCAL_BEAM(:,2),2)).^2+(COORD(LOCAL_BEAM(:,3),3)-COORD(LOCAL_BEAM(:,2),3)).^2));
A=1;
I=1;
J=2;
nu=0.3;
G=E./2./(1+nu);
for k=1:size(Beam_Element_connectivity,1)
    c1=cos(theta(k));
    s1=sin(theta(k));
    c2=cos(-psi(k));
    s2=sin(-psi(k));
    R1=[c1 s1 0
        -s1 c1 0
        0 0 1];
    R2=[c2 0 -s2
        0 1 0
        s2 0 c2];
    T=R2*R1;
    Z=zeros(3,3);
    T=[ T Z Z Z
        Z T Z Z
        Z Z T Z
        Z Z Z T];
    %local stiffness matrix
    a=Beam_Length(k)/2;
    Kt=A*E(k)/a/2;
    Kf=E(k)*I/a^3/2;
    Kr=G(k)*J/a/2;
    m=12;
    line_vector=[1;1;2;2;2;2;3;3;3;3;4;4;5;5;5;6;6;6;7;8;8;9;9;10;11;12];
    column_vector=[1;7;2;6;8;12;3;5;9;11;4;10;5;9;11;6;8;12;7;8;12;9;11;10;11;12];
    Component=[Kt;-Kt;3*Kf;3*a*Kf;-3*Kf;3*a*Kf;3*Kf;-3*a*Kf;-3*Kf;-3*a*Kf;Kr;-Kr;4*a^2*Kf;3*Kf*a;Kf*a^2*2;4*a^2*Kf;-3*a*Kf;2*a^2*Kf;Kt;3*Kf;-3*a*Kf;3*Kf;3*a*Kf;Kr;4*a^2*Kf;4*a^2*Kf];
    Kel_local=sparse(line_vector,column_vector,Component,m,m);
    Kel_local=Kel_local+tril(Kel_local.',-1);
    KEl_local(k).mat=Kel_local;
    Kel_global(:,:,k)=T.'*Kel_local*T;
end
%truss matrix assembly
At=3066;
theta_rod=atan2(COORD(LOCAL_ROD(:,3),3)-COORD(LOCAL_ROD(:,2),3),COORD(LOCAL_ROD(:,3),2)-COORD(LOCAL_ROD(:,2),2));
psi_rod=atan2(COORD(LOCAL_ROD(:,3),4)-COORD(LOCAL_ROD(:,2),4),sqrt((COORD(LOCAL_ROD(:,3),2)-COORD(LOCAL_ROD(:,2),2)).^2+(COORD(LOCAL_ROD(:,3),3)-COORD(LOCAL_ROD(:,2),3)).^2));
for k=1:size(Rod_connectivity_matrix,1)
    c1=cos(theta_rod(k));
    s1=sin(theta_rod(k));
    c2=cos(-psi_rod(k));
    s2=sin(-psi_rod(k));
    R1=[c1 s1 0
        -s1 c1 0
        0 0 1];
    R2=[c2 0 -s2
        0 1 0
        s2 0 c2];
    T=R2*R1;
    Z=zeros(3,3);
    T=[ T Z Z Z
        Z T Z Z
        Z Z T Z
        Z Z Z T];
    %local stiffness matrix
    a=Rod_Length(k)/2;
    Kt=At*Erod/a/2;
    Kf=0;
    Kr=0;
    m=12;
    line_vector=[1;1;2;2;2;2;3;3;3;3;4;4;5;5;5;6;6;6;7;8;8;9;9;10;11;12];
    column_vector=[1;7;2;6;8;12;3;5;9;11;4;10;5;9;11;6;8;12;7;8;12;9;11;10;11;12];
    Component=[Kt;-Kt;3*Kf;3*a*Kf;-3*Kf;3*a*Kf;3*Kf;-3*a*Kf;-3*Kf;-3*a*Kf;Kr;-Kr;4*a^2*Kf;3*Kf*a;Kf*a^2*2;4*a^2*Kf;-3*a*Kf;2*a^2*Kf;Kt;3*Kf;-3*a*Kf;3*Kf;3*a*Kf;Kr;4*a^2*Kf;4*a^2*Kf];
    Kel_local=sparse(line_vector,column_vector,Component,m,m);
    Kel_local=Kel_local+tril(Kel_local.',-1);
    KEl_local_ROD(k).mat=Kel_local;
    Kel_global_ROD(:,:,k)=T.'*Kel_local*T;
end
Number_of_DOF=max(DOFs_MAP(:,4));
K=zeros(Number_of_DOF);
F=zeros(Number_of_DOF,3);
% K(1:Number_of_retained_DOF,1:Number_of_retained_DOF)=K_interface;
K(1:Number_of_retained_DOF,1:Number_of_retained_DOF)=K_interface;
F(1:Number_of_retained_DOF,:)=LOAD_MATRIX;
%Stiffness assembly (once for each configuration)
for k=1:size(Beam_Element_connectivity,1)+size(Rod_connectivity_matrix,1)
    if k<=size(Beam_Element_connectivity,1)
        K=K+recontstruction_matrix(:,:,k)*Kel_global(:,:,k)*recontstruction_matrix(:,:,k).';
    else
        K=K+recontstruction_matrix(:,:,k)*Kel_global_ROD(:,:,k-size(Beam_Element_connectivity,1))*recontstruction_matrix(:,:,k).';
    end
end
%BCs application
%DOFs 123456 in node P2 and P2d
K=K(DOFs_MAP(:,1)~=P2(1)&DOFs_MAP(:,1)~=P2d(1),DOFs_MAP(:,1)~=P2(1)&DOFs_MAP(:,1)~=P2d(1));
disp('Stiffness Matrix Assembly Ended')
toc