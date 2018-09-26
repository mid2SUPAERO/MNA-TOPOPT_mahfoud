%connectivity matrix
clear
close
load('lin_pert_sub.dat.mat')
x_initial=4.406290039E+03;
x=4500;
Radius=norm([668.671021,121.383797]);
initial_phase=atan2(668.671021,121.383797);
% T=initial_phase*180/pi;
T=20;
P1=[(max(intercoord(:,1)+1)) 2.451298096E+03,-Radius*sind(T),Radius*cosd(T)];
P2=[(max(intercoord(:,1)+2)) x,-66.380302,910.036987];
P1d=[(max(intercoord(:,1)+3)) 2.451298096E+03,Radius*sind(T),Radius*cosd(T)];
P2d=[(max(intercoord(:,1)+4)) x,66.380302,910.036987];
P=[P1;P2;P1d;P2d];
right_nodes=intercoord(intercoord(:,3)>0,:);
left_nodes=intercoord(intercoord(:,3)<0,:);
Beam_Element_connectivity=zeros((size(intercoord,1)),3);
Beam_Element_connectivity(:,1)=(1:size(intercoord,1)).';
Beam_Element_connectivity(:,2)=[P(1,1)*ones(length(left_nodes),1);P(3,1)*ones(length(right_nodes),1)];
Beam_Element_connectivity(:,3)=[left_nodes(:,1);right_nodes(:,1)];
Rod_connectivity_matrix=[size(intercoord,1)+1 P1(1) P2(1);size(intercoord,1)+2 P1d(1) P2d(1)];
COORD=[intercoord;P];
load('connectivity_interface');
LOCAL_QUAD=ELint;
for k=2:5
    for l=1:size(ELint,1)
        LOCAL_QUAD(l,k)=find(COORD(:,1)==ELint(l,k));
    end
end
LOCAL_BEAM=Beam_Element_connectivity;
for k=2:3
    for l=1:size(Beam_Element_connectivity,1)
        LOCAL_BEAM(l,k)=find(COORD(:,1)==Beam_Element_connectivity(l,k));
    end
end
LOCAL_ROD=Rod_connectivity_matrix;
for k=2:3
    for l=1:size(Rod_connectivity_matrix,1)
        LOCAL_ROD(l,k)=find(COORD(:,1)==Rod_connectivity_matrix(l,k));
    end
end
%Calculation of Beam and Rod length
Beam_Length=sqrt((COORD(LOCAL_BEAM(:,3),2)-COORD(LOCAL_BEAM(:,2),2)).^2+(COORD(LOCAL_BEAM(:,3),3)-COORD(LOCAL_BEAM(:,2),3)).^2+(COORD(LOCAL_BEAM(:,3),4)-COORD(LOCAL_BEAM(:,2),4)).^2);
Rod_Length=sqrt((COORD(LOCAL_ROD(:,3),2)-COORD(LOCAL_ROD(:,2),2)).^2+(COORD(LOCAL_ROD(:,3),3)-COORD(LOCAL_ROD(:,2),3)).^2+(COORD(LOCAL_ROD(:,3),4)-COORD(LOCAL_ROD(:,2),4)).^2);
phase_square_difference=(atan2(COORD(LOCAL_BEAM(:,2),4),COORD(LOCAL_BEAM(:,2),3))-atan2(COORD(LOCAL_BEAM(:,3),4),COORD(LOCAL_BEAM(:,3),3))).^2;
%Approximation of Heaviside function for Young module evaluation
betha=50;threshold=0.05;
H=@(x) 1-1./(1+exp(-2*betha*(x-threshold)));
%Young Module Evaluation
Emax=1e9;
E=Emax*H(phase_square_difference);
%Rod Young Module
Erod=210000;
Display = 1;
if (Display == 1)
    load('lin_pert_sub.dat.mat')
    
    h = figure(1); set(h,'Color',[1 1 1]);
    
    hold on; h24_patch = patch('Vertices',COORD(:,2:4),'Faces',LOCAL_QUAD(:,2:5),'FaceVertexCData',ones(size(LOCAL_QUAD,1),1)*[0.3 0.1 0.4],'FaceColor','flat');
    for k=1:size(Rod_connectivity_matrix,1)
        plot3(COORD(LOCAL_ROD(k,2:3),2),COORD(LOCAL_ROD(k,2:3),3),COORD(LOCAL_ROD(k,2:3),4),'b','LineWidth',3)
    end
    for k=1:size(Beam_Element_connectivity,1)
        plot3(COORD(LOCAL_BEAM(k,2:3),2),COORD(LOCAL_BEAM(k,2:3),3),COORD(LOCAL_BEAM(k,2:3),4),'Color',[H(phase_square_difference(k)),1-H(phase_square_difference(k)),0],'LineWidth',5*(H(phase_square_difference(k))*(1-0.01)+0.01))
    end
    view([51.52 24.72]);
    axis equal;
    axis off;
    figure(2)
    plot(2*betha*phase_square_difference,E,'o',2*betha*(0:0.001:max(phase_square_difference)),Emax*H((0:0.001:max(phase_square_difference))),'LineWidth',3)
    grid on
    xlabel('x = 2*\beta*(\theta - \theta_{c})^2')
    ylabel('Beam Young Module')
    legend('Beam Young Module','E(x) = 1-1/(1+exp(-(x-x_{th})))')
end
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
load('substructure_results.mat')
Number_of_retained_DOF=max(DOFs_MAP(:,4));
for k=1:size(P,1)
    DOFs_MAP=[DOFs_MAP;[P(k,1)*ones(6,1),(max(DOFs_MAP(:,2))+1)*ones(6,1),(1:6).',(max(DOFs_MAP(:,4)))+(1:6).']];
end
Number_of_DOF=max(DOFs_MAP(:,4));
K=zeros(Number_of_DOF);
F=zeros(Number_of_DOF,3);
K(1:Number_of_retained_DOF,1:Number_of_retained_DOF)=K_interface;
F(1:Number_of_retained_DOF,:)=LOAD_MATRIX;
%reconstruction matrix(to be done once for all)
recontstruction_matrix=zeros(Number_of_DOF,12,size(Beam_Element_connectivity,1)+size(Rod_connectivity_matrix,1));
for k=1:size(Beam_Element_connectivity,1)+size(Rod_connectivity_matrix,1)
    for n=2:3
        for dof=1:6
            if k<=size(Beam_Element_connectivity,1)
                recontstruction_matrix((DOFs_MAP(:,1)==Beam_Element_connectivity(k,n))&(DOFs_MAP(:,3)==dof),6*(n-2)+dof,k)=1;
            else
                recontstruction_matrix((DOFs_MAP(:,1)==Rod_connectivity_matrix(k-size(Beam_Element_connectivity,1),n))&(DOFs_MAP(:,3)==dof),6*(n-2)+dof,k)=1;
            end
        end
    end
end
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
F=F(DOFs_MAP(:,1)~=P2(1)&DOFs_MAP(:,1)~=P2d(1),:);
Final_DOF_MAP=DOFs_MAP(DOFs_MAP(:,1)~=P2(1)&DOFs_MAP(:,1)~=P2d(1),:);
%Displacement evaluation
tic
U=K\F;
toc
%Disp-tip clearance
U_out=U(1:Number_of_retained_DOF,:);
Utip=Recovery_Matrix*U_out;
fichier='connectivity_mat';
eval([ 'save ' fichier ' Utip recontstruction_matrix Beam_Element_connectivity Rod_connectivity_matrix LOCAL_BEAM LOCAL_ROD DOFs_MAP Final_DOF_MAP Number_of_retained_DOF']);
disp([ 'LE FICHIER ' fichier '.mat A ETE CREE']);
if (Display == 1)
    New_COORD=COORD;
    scalefact=1;
    for k=1:size(COORD,1)
        for dof=1:3
            if any(Final_DOF_MAP(:,1)==New_COORD(k,1))
                New_COORD(k,dof+1)=New_COORD(k,dof+1)+scalefact*U(Final_DOF_MAP(:,1)==New_COORD(k,1)&Final_DOF_MAP(:,3)==dof,1);
            end
        end
    end
    
    h = figure(1); set(h,'Color',[1 1 1]);
    
    hold on; h24_patch = patch('Vertices',New_COORD(:,2:4),'Faces',LOCAL_QUAD(:,2:5),'FaceVertexCData',ones(size(LOCAL_QUAD,1),1)*[0.3 0.1 0.4],'FaceColor','flat');
    for k=1:size(Rod_connectivity_matrix,1)
        plot3(New_COORD(LOCAL_ROD(k,2:3),2),New_COORD(LOCAL_ROD(k,2:3),3),New_COORD(LOCAL_ROD(k,2:3),4),'b','LineWidth',3)
    end
    for k=1:size(Beam_Element_connectivity,1)
        plot3(New_COORD(LOCAL_BEAM(k,2:3),2),New_COORD(LOCAL_BEAM(k,2:3),3),New_COORD(LOCAL_BEAM(k,2:3),4),'Color',[H(phase_square_difference(k)),1-H(phase_square_difference(k)),0],'LineWidth',5*(H(phase_square_difference(k))*(1-0.01)+0.01))
    end
    view([51.52 24.72]);
    axis equal;
    axis off;
    figure(2)
    plot(2*betha*phase_square_difference,E,'o',2*betha*(0:0.001:max(phase_square_difference)),Emax*H((0:0.001:max(phase_square_difference))),'LineWidth',3)
    grid on
    xlabel('x = 2*\beta*(\theta - \theta_{c})^2')
    ylabel('Beam Young Module')
    legend('Beam Young Module','E(x) = 1-1/(1+exp(-(x-x_{th})))')
end



