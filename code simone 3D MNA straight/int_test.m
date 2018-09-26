%% TEST INTERNODES
%%
% *Geometry definition*
%
%
% Cube geometry:
%
% $$ \left\lbrace \begin{array}{c} -1\leq x \leq 1 \\ -1 \leq y \leq 1 \\ 0 \leq z \leq 2 \end{array}\right. $$
%
% Plate geometry:
%
% $$ \left\lbrace \begin{array}{c} -2\leq x \leq 2 \\ -2 \leq y \leq 2 \\ -1 \leq z \leq 0 \end{array}\right. $$
%
% *mesh definition*
%
% $$ n_x^c $$ : number of elements on x direcrion
%
% $$ n_y^c $$ : number of elements on y direcrion
%
% $$ n_z^c $$ : number of elements on z direcrion
clear
close all
% clear
% Nn=4;
% nxc=3*Nn;
% nyc=3*Nn;
% nzc=3*Nn;
% nxp=2*Nn;
% nyp=2*Nn;
% nzp=2*Nn;
% lxc=2;
% lyc=2;
% lzc=2;
% lxp=2;
% lyp=2;
% lzp=2;
conf=1;
lc=2;
correzione_momento=2;
Intergrid_switch=16;
if lc==1
Load_dir=3;
else
Load_dir=[1,3];
end
if conf==1
%%conf 1
nxc=10;
nyc=10;
nzc=10;
nxp=4;
nyp=nxp;
nzp=20;
lxc=2;
lyc=2;
lzc=2;
lxp=2;
lyp=2;
lzp=4;
else
%%conf 2
nxc=10;
nyc=10;
nzc=10;
nxp=20;
nyp=nxp;
nzp=5;
lxc=2;
lyc=2;
lzc=2;
lxp=4;
lyp=4;
lzp=1;
end
interface_cube=2;

apriori_patch=0;
pressure_switch=1;
primal_schwartz=0;
Vhon_mises=1;
nomaster=0;
cube_master=1;
RL_RBF=1;

curbature=3;
average_stress=1; global_least_sq=1;
interp_residual=0;
xc=linspace(-lxc/2,lxc/2,nxc+1);
yc=linspace(-lyc/2,lyc/2,nyc+1);
zc=linspace(0,lzc,nzc+1);
xp=linspace(-lxp/2,lxp/2,nxp+1);
yp=linspace(-lyp/2,lyp/2,nyp+1);
zp=linspace(-lzp,0,nzp+1);
[Xc,Yc,Zc]=meshgrid(xc,yc,zc);
[Xp,Yp,Zp]=meshgrid(xp,yp,zp);
Xc=Xc(:);Xp=Xp(:);
Yc=Yc(:);Yp=Yp(:);
Zc=Zc(:);Zp=Zp(:);
interface_c=find(Zc==0);
not_int_c=setdiff(1:length(Zc),interface_c);
if interface_cube==1
    interface_p=find(Zp==0&abs(Xp)<=lxc/2&abs(Yp)<=lyc/2);
elseif interface_cube==0
    interface_p=find(Zp==0);
    
else
    interface_p=find(Zp==0);
    internal_nodes=find(Zp==0&abs(Xp)<lxc/2&abs(Yp)<lyc/2);
end

external_node=find(Zp==0&(abs(Xp)>lxc/2|abs(Yp)>lyc/2));
[~,external_node_id]=intersect(interface_p,external_node);
clamped_nodes=find(Zp==min(Zp));
total_clamped=find(Zp==min(Zp)&Xp==min(Xp)&Yp==min(Yp));
yzclamped=find(Zp==min(Zp)&Xp==max(Xp)&Yp==min(Yp));
excited_nodes=find(Zc==lzc);
if curbature==1;
    Z0=-1.5*lzp;
    rc=Zc-Z0;
    rp=Zp-Z0;
    tethac=Xc./rc;
    tethap=Xp./rp;
    Xp=rp.*sin(tethap);
    Zp=Z0+rp.*(cos(tethap));
    Xc=rc.*sin(tethac);
    Zc=Z0+rc.*(cos(tethac));
    %centre constante
elseif curbature==2;
    %curbure constante
    Z0=-0.5*lzp;
    rc=-Z0;
    rp=-Z0;
    tethac=Xc./rc;
    tethap=Xp./rp;
    Xp=rp.*sin(tethap);
    Zp=Zp+Z0+rp.*(cos(tethap));
    Xc=rc.*sin(tethac);
    Zc=Zc+Z0+rc.*(cos(tethac));
elseif curbature==3;
    %double curbure constante
    Z0=-0.5*4;
    rc=-Z0;
    rp=-Z0;
    tethac=Xc./rc;
    phic=Yc./rc;
    tethap=Xp./rp;
    phip=Yp./rp;
    Zp=Zp+Z0+rp./(sqrt(1+sin(phip).^2+sin(tethap).^2));
    Xp=rp./(sqrt(1+sin(phip).^2+sin(tethap).^2)).*sin(tethap);
    Yp=rp./(sqrt(1+sin(phip).^2+sin(tethap).^2)).*sin(phip);
    Zc=Zc+Z0+rc./(sqrt(1+sin(phic).^2+sin(tethac).^2));
    Xc=rc./(sqrt(1+sin(phic).^2+sin(tethac).^2)).*sin(tethac);
    Yc=rc./(sqrt(1+sin(phic).^2+sin(tethac).^2)).*sin(phic);
    %     Xp=rp.*sin(tethap);
    %     Zp=Zp+Z0+rp.*(cos(tethap));
    %     Xc=rc.*sin(tethac);
    %     Zc=Zc+Z0+rc.*(cos(tethac));
elseif curbature==4;
    xm=lxp/2;ym=lyp/2;
    DZm=0.5;
    Zp=Zp+DZm*(1-Zp.^2/(lzp*lzc)+(lzc-lzp)/(lzp*lzc)*Zp)/(xm^2+ym^2).*(xm^2-Xp.^2+ym.^2-Yp.^2);
    Zc=Zc+DZm*(1-Zc.^2/(lzp*lzc)+(lzc-lzp)/(lzp*lzc)*Zc)/(xm^2+ym^2).*(xm^2-Xc.^2+ym.^2-Yc.^2);
elseif curbature==5;
    xm=lxp/2;ym=lyp/2;
    DZm=0.5;
    Zp=Zp+DZm*(1-Zp.^2/(lzp*2*lzc)+(2*lzc-lzp)/(2*lzp*lzc)*Zp)/(xm^2+ym^2).*(xm^2-Xp.^2).*(ym.^2-Yp.^2);
    Zc=Zc+DZm*(1-Zc.^2/(lzp*2*lzc)+(2*lzc-lzp)/(lzp*2*lzc)*Zc)/(xm^2+ym^2).*(xm^2-Xc.^2).*(ym.^2-Yc.^2);
    
end
Cube_coord=[Xc,Yc,Zc];
Plate_coord=[Xp,Yp,Zp];
figure(1)
scatter3([Xc;Xp],[Yc;Yp],[Zc;Zp],'fil','r')
axis equal
title('Cube Grid')
xlabel x
ylabel y
zlabel z
Nc=length(Xc); Np=length(Xp);
% for m=1:Nc
%     text(Xc(m),Yc(m),Zc(m),num2str(m))
% end
% for m=1:Np
%     text(Xp(m),Yp(m),Zp(m),num2str(m+Nc))
% end
%Element connectivity
Nec=nxc*nyc*nzc; Nep=nxp*nyp*nzp;
Ec=zeros(Nec,8);
for n=1:Nec
    Ec(n,:)=[n n+nyc+1 n+nyc+2 n+1 n+(nyc+1)*(nxc+1) n+(nyc+1)*(nxc+1)+nyc+1 n+(nyc+1)*(nxc+1)+nyc+2 n+(nyc+1)*(nxc+1)+1]+fix((n-1)/(nyc))+(nyc+1)*fix((n-1)/(nxc*nyc));
end
%Faces
FACES_c=zeros(6*Nec,5);
FACES_c(:,1)=1:6*Nec;
for k=1:Nec
    FACES_c(6*(k-1)+1,2:end)=Ec(k,[1 2 3 4]);
    FACES_c(6*(k-1)+2,2:end)=Ec(k,[1 2 6 5]);
    FACES_c(6*(k-1)+3,2:end)=Ec(k,[2 3 7 6]);
    FACES_c(6*(k-1)+4,2:end)=Ec(k,[5 6 7 8]);
    FACES_c(6*(k-1)+5,2:end)=Ec(k,[1 4 8 5]);
    FACES_c(6*(k-1)+6,2:end)=Ec(k,[4 3 7 8]);
end
Ep=zeros(Nep,8);
for n=1:Nep
    Ep(n,:)=[n n+nyp+1 n+nyp+2 n+1 n+(nyp+1)*(nxp+1)  n+(nyp+1)*(nxp+1)+nyp+1 n+(nyp+1)*(nxp+1)+nyp+2 n+(nyp+1)*(nxp+1)+1]+fix((n-1)/(nyp))+(nyp+1)*fix((n-1)/(nxp*nyp));
end
%Faces
FACES_p=zeros(6*Nep,5);
FACES_p(:,1)=1:6*Nep;
for k=1:Nep
    FACES_p(6*(k-1)+1,2:end)=Ep(k,[1 2 3 4]);
    FACES_p(6*(k-1)+2,2:end)=Ep(k,[1 2 6 5]);
    FACES_p(6*(k-1)+3,2:end)=Ep(k,[2 3 7 6]);
    FACES_p(6*(k-1)+4,2:end)=Ep(k,[5 6 7 8]);
    FACES_p(6*(k-1)+5,2:end)=Ep(k,[1 4 8 5]);
    FACES_p(6*(k-1)+6,2:end)=Ep(k,[4 3 7 8]);
end
%display the mesh
figure(2)
hold on; h24_patch = patch('Vertices',Cube_coord,'Faces',FACES_c(:,2:5),'FaceVertexCData',[0 1 1],'FaceColor','flat','LineWidth',1); axis equal; axis off; drawnow;
hold on; h24_patch = patch('Vertices',Plate_coord,'Faces',FACES_p(:,2:5),'FaceVertexCData',[0 0 1],'FaceColor','flat','LineWidth',1); axis equal; axis off; drawnow;
quiver3(0,0,0,4,0,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,4,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,0,4,'LineWidth',3,'Color','k')
text(4,0,0,'x','FontSize',18,'FontWeight','bold')
text(0,4,0,'y','FontSize',18,'FontWeight','bold')
text(0,0,4,'z','FontSize',18,'FontWeight','bold')
view([135 35]);
hold off
% dofs identification cube
first_node_index_c=Ec(:,1);
second_node_index_c=Ec(:,2);
third_node_index_c=Ec(:,3);
fourth_node_index_c=Ec(:,4);
fivth_node_index_c=Ec(:,5);
sixth_node_index_c=Ec(:,6);
seventh_node_index_c=Ec(:,7);
eigth_node_index_c=Ec(:,8);
first_node_DOFs_index_c=[3*(first_node_index_c-1)+1,3*(first_node_index_c-1)+2,3*(first_node_index_c-1)+3];
second_node_DOFs_index_c=[3*(second_node_index_c-1)+1,3*(second_node_index_c-1)+2,3*(second_node_index_c-1)+3];
third_node_DOFs_index_c=[3*(third_node_index_c-1)+1,3*(third_node_index_c-1)+2,3*(third_node_index_c-1)+3];
fourth_node_DOFs_index_c=[3*(fourth_node_index_c-1)+1,3*(fourth_node_index_c-1)+2,3*(fourth_node_index_c-1)+3];
fivth_node_DOFs_index_c=[3*(fivth_node_index_c-1)+1,3*(fivth_node_index_c-1)+2,3*(fivth_node_index_c-1)+3];
sixth_node_DOFs_index_c=[3*(sixth_node_index_c-1)+1,3*(sixth_node_index_c-1)+2,3*(sixth_node_index_c-1)+3];
seventh_node_DOFs_index_c=[3*(seventh_node_index_c-1)+1,3*(seventh_node_index_c-1)+2,3*(seventh_node_index_c-1)+3];
eighth_node_DOFs_index_c=[3*(eigth_node_index_c-1)+1,3*(eigth_node_index_c-1)+2,3*(eigth_node_index_c-1)+3];
CUBE_ZONE_ELEMENT_DOFS=[first_node_DOFs_index_c,second_node_DOFs_index_c,third_node_DOFs_index_c,fourth_node_DOFs_index_c,fivth_node_DOFs_index_c,sixth_node_DOFs_index_c,seventh_node_DOFs_index_c,eighth_node_DOFs_index_c];
Cube_DOFs_Map=zeros(3*size(Cube_coord,1),4);
for k=1:size(Cube_coord,1)
    Cube_DOFs_Map(3*k-2,:)=[k,k,1,3*k-2];
    Cube_DOFs_Map(3*k-1,:)=[k,k,2,3*k-1];
    Cube_DOFs_Map(3*k,:)=[k,k,3,3*k];
end
% dofs identification plate
first_node_index_p=Ep(:,1);
second_node_index_p=Ep(:,2);
third_node_index_p=Ep(:,3);
fourth_node_index_p=Ep(:,4);
fivth_node_index_p=Ep(:,5);
sixth_node_index_p=Ep(:,6);
seventh_node_index_p=Ep(:,7);
eigth_node_index_p=Ep(:,8);
first_node_DOFs_index_p=[3*(first_node_index_p-1)+1,3*(first_node_index_p-1)+2,3*(first_node_index_p-1)+3];
second_node_DOFs_index_p=[3*(second_node_index_p-1)+1,3*(second_node_index_p-1)+2,3*(second_node_index_p-1)+3];
third_node_DOFs_index_p=[3*(third_node_index_p-1)+1,3*(third_node_index_p-1)+2,3*(third_node_index_p-1)+3];
fourth_node_DOFs_index_p=[3*(fourth_node_index_p-1)+1,3*(fourth_node_index_p-1)+2,3*(fourth_node_index_p-1)+3];
fivth_node_DOFs_index_p=[3*(fivth_node_index_p-1)+1,3*(fivth_node_index_p-1)+2,3*(fivth_node_index_p-1)+3];
sixth_node_DOFs_index_p=[3*(sixth_node_index_p-1)+1,3*(sixth_node_index_p-1)+2,3*(sixth_node_index_p-1)+3];
seventh_node_DOFs_index_p=[3*(seventh_node_index_p-1)+1,3*(seventh_node_index_p-1)+2,3*(seventh_node_index_p-1)+3];
eighth_node_DOFs_index_p=[3*(eigth_node_index_p-1)+1,3*(eigth_node_index_p-1)+2,3*(eigth_node_index_p-1)+3];
PLATE_ZONE_ELEMENT_DOFS=[first_node_DOFs_index_p,second_node_DOFs_index_p,third_node_DOFs_index_p,fourth_node_DOFs_index_p,fivth_node_DOFs_index_p,sixth_node_DOFs_index_p,seventh_node_DOFs_index_p,eighth_node_DOFs_index_p];
Plate_DOFs_Map=zeros(3*size(Plate_coord,1),4);
for k=1:size(Plate_coord,1)
    Plate_DOFs_Map(3*k-2,:)=[k,k,1,3*k-2];
    Plate_DOFs_Map(3*k-1,:)=[k,k,2,3*k-1];
    Plate_DOFs_Map(3*k,:)=[k,k,3,3*k];
end
%GLOBAL DOFs
Global_DOFs_MAP=[[Cube_DOFs_Map(:,1);Plate_DOFs_Map(:,1)+Cube_DOFs_Map(end,1)],[Cube_DOFs_Map(:,2:3);Plate_DOFs_Map(:,2:3)],[Cube_DOFs_Map(:,4);Plate_DOFs_Map(:,4)+Cube_DOFs_Map(end,4)]];
%cube interface mass matrix

int_dofs_c=zeros(3*length(interface_c),5);
for k=1:size(interface_c,1)
    int_dofs_c(3*k-2,:)=[k,interface_c(k),1,3*k-2,Cube_DOFs_Map(interface_c(k)==Cube_DOFs_Map(:,1)&(1==Cube_DOFs_Map(:,3)),4)];
    int_dofs_c(3*k-1,:)=[k,interface_c(k),2,3*k-1,Cube_DOFs_Map(interface_c(k)==Cube_DOFs_Map(:,1)&(2==Cube_DOFs_Map(:,3)),4)];
    int_dofs_c(3*k,:)=[k,interface_c(k),3,3*k,Cube_DOFs_Map(interface_c(k)==Cube_DOFs_Map(:,1)&(3==Cube_DOFs_Map(:,3)),4)];
end
interface_c_coord=Cube_coord(interface_c,:);
[ID4]=ismember(FACES_c(:,[2:5]),interface_c);
IDF=ID4(:,1)&ID4(:,2)&ID4(:,3)&ID4(:,4);
FACE_int_c=FACES_c(IDF,:);
[Mc]=interface_mass_assembly(FACE_int_c,Cube_DOFs_Map,Cube_coord);
Mc=Mc(int_dofs_c(:,5),int_dofs_c(:,5));
%plate interface mass matrix


int_dofs_p=zeros(3*length(interface_p),5);
for k=1:size(interface_p,1)
    int_dofs_p(3*k-2,:)=[k,interface_p(k),1,3*k-2,Plate_DOFs_Map(interface_p(k)==Plate_DOFs_Map(:,1)&(1==Plate_DOFs_Map(:,3)),4)];
    int_dofs_p(3*k-1,:)=[k,interface_p(k),2,3*k-1,Plate_DOFs_Map(interface_p(k)==Plate_DOFs_Map(:,1)&(2==Plate_DOFs_Map(:,3)),4)];
    int_dofs_p(3*k,:)=[k,interface_p(k),3,3*k,Plate_DOFs_Map(interface_p(k)==Plate_DOFs_Map(:,1)&(3==Plate_DOFs_Map(:,3)),4)];
end
interface_p_coord=Plate_coord(interface_p,:);
[ID4]=ismember(FACES_p(:,[2:5]),interface_p);
IDF=ID4(:,1)&ID4(:,2)&ID4(:,3)&ID4(:,4);
FACE_int_p=FACES_p(IDF,:);
if interface_cube==2
    [ID4]=ismember(FACE_int_p(:,[2:5]),internal_nodes);
    IDF=ID4(:,1)|ID4(:,2)|ID4(:,3)|ID4(:,4);
    FACE_int_p=FACE_int_p(IDF,:);
    interface_p=unique(reshape(FACE_int_p(:,2:end),[],1));
    not_int_p=setdiff(1:length(Zp),interface_p);
    interface_p_coord=Plate_coord(interface_p,:);
    int_dofs_p=zeros(3*length(interface_p),5);
    for k=1:size(interface_p,1)
        int_dofs_p(3*k-2,:)=[k,interface_p(k),1,3*k-2,Plate_DOFs_Map(interface_p(k)==Plate_DOFs_Map(:,1)&(1==Plate_DOFs_Map(:,3)),4)];
        int_dofs_p(3*k-1,:)=[k,interface_p(k),2,3*k-1,Plate_DOFs_Map(interface_p(k)==Plate_DOFs_Map(:,1)&(2==Plate_DOFs_Map(:,3)),4)];
        int_dofs_p(3*k,:)=[k,interface_p(k),3,3*k,Plate_DOFs_Map(interface_p(k)==Plate_DOFs_Map(:,1)&(3==Plate_DOFs_Map(:,3)),4)];
    end
    [~,external_node_id]=intersect(interface_p,external_node);
end
[Mp]=interface_mass_assembly(FACE_int_p,Plate_DOFs_Map,Plate_coord);
Mp=Mp(int_dofs_p(:,5),int_dofs_p(:,5));
%Projection operators RL-RBF

if RL_RBF==1
    %DZ->E
    DISTANCE_cp=zeros(size(interface_c_coord,1),size(interface_p_coord,1));
    for n=1:size(interface_p_coord,1)
        DISTANCE_cp(:,n)=sqrt(diag((interface_c_coord-ones(size(interface_c_coord,1),1)*interface_p_coord(n,:))*(interface_c_coord-ones(size(interface_c_coord,1),1)*interface_p_coord(n,:)).'));
    end
    DISTANCE_cc=zeros(size(interface_c_coord,1));
    for n=1:size(interface_c_coord,1)
        DISTANCE_cc(:,n)=sqrt(diag((interface_c_coord-ones(size(interface_c_coord,1),1)*interface_c_coord(n,:))*(interface_c_coord-ones(size(interface_c_coord,1),1)*interface_c_coord(n,:)).'));
    end
    DISTANCE_pp=zeros(size(interface_p_coord,1),size(interface_p_coord,1));
    for n=1:size(interface_p_coord,1)
        DISTANCE_pp(:,n)=sqrt(diag((interface_p_coord-ones(size(interface_p_coord,1),1)*interface_p_coord(n,:))*(interface_p_coord-ones(size(interface_p_coord,1),1)*interface_p_coord(n,:)).'));
    end
    
    %reference Radii evaluation
    connectivity_matrix_c=zeros(9,size(interface_c_coord,1));
    Neig_dist_c=zeros(9,size(interface_c_coord,1));
    for nod=1:size(interface_c_coord,1)
        [Iel,~]=find(nod==FACE_int_c(:,2:5));
        [neighbourhood_el]=(FACE_int_c(Iel,2:5));
        neighbourhood_el=unique(neighbourhood_el(:));
        neighbourhood_el=setdiff(neighbourhood_el,nod);
        connectivity_matrix_c(1:(length(neighbourhood_el)+1),nod)=[nod;neighbourhood_el];
        for l=1:9
            if connectivity_matrix_c(l,nod)~=0
                Neig_dist_c(l,nod)=DISTANCE_cc(nod,connectivity_matrix_c(l,nod));
            else
                Neig_dist_c(l,nod)=0;
            end
        end
    end
    
    ELintip=FACE_int_p;
    for k=1:size(FACE_int_p,1)
        for l=1:size(FACE_int_p,2)-1
            ELintip(k,l+1)=find(interface_p==FACE_int_p(k,l+1));
        end
    end
    ELintic=FACE_int_c;
    for k=1:size(FACE_int_c,1)
        for l=1:size(FACE_int_c,2)-1
            ELintic(k,l+1)=find(interface_c==FACE_int_c(k,l+1));
        end
    end
    connectivity_matrix_p=zeros(9,size(interface_p_coord,1));
    Neig_dist_p=zeros(9,size(interface_p_coord,1));
    for nod=1:size(interface_p_coord,1)
        [Iel,~]=find((nod)==ELintip(:,2:5));
        [neighbourhood_el]=(ELintip(Iel,2:5));
        neighbourhood_el=unique(neighbourhood_el(:));
        neighbourhood_el=setdiff(neighbourhood_el,(nod));
        connectivity_matrix_p(1:(length(neighbourhood_el)+1),nod)=[(nod);neighbourhood_el];
        for l=1:9
            if connectivity_matrix_p(l,nod)~=0
                Neig_dist_p(l,nod)=DISTANCE_pp(nod,connectivity_matrix_p(l,nod));
            else
                Neig_dist_p(l,nod)=0;
            end
        end
    end
    radii_p=max(Neig_dist_p);
    radii_c=max(Neig_dist_c);
    radii_p=repmat(max(max(radii_p),max(radii_c)),size(radii_p,1),size(radii_p,2));
    radii_c=repmat(max(max(radii_p),max(radii_c)),size(radii_c,1),size(radii_c,2));
    radii_pp=repmat(radii_p,size(DISTANCE_pp,1),1);
    PHI_MM_p=(1-DISTANCE_pp./radii_pp).^4.*(1-DISTANCE_pp./radii_pp>=0).*(1+4.*DISTANCE_pp./radii_pp);
    radii_cp=repmat(radii_p,size(DISTANCE_cc,1),1);
    PHI_NM_p=(1-DISTANCE_cp./radii_cp).^4.*(1-DISTANCE_cp./radii_cp>=0).*(1+4.*DISTANCE_cp./radii_cp);
    radii_cc=repmat(radii_c,size(DISTANCE_cc,1),1);
    if interface_cube==0||interface_cube==2
        % pressure of external nodes is anyway 0
        PHI_NM_p_pressure=PHI_NM_p;
        PHI_NM_p_pressure(:,external_node_id)=0;
    end
    PHI_MM_c=(1-DISTANCE_cc./radii_cc).^4.*(1-DISTANCE_cc./radii_cc>=0).*(1+4.*DISTANCE_cc./radii_cc);
    radii_pc=repmat(radii_c,size(DISTANCE_cp,2),1);
    PHI_NM_c=(1-DISTANCE_cp.'./radii_pc).^4.*(1-DISTANCE_cp.'./radii_pc>=0).*(1+4.*DISTANCE_cp.'./radii_pc);
    if interface_cube==0||interface_cube==2
        % pressure of external nodes is anyway 0
        PHI_NM_c_pressure=PHI_NM_c;
        PHI_NM_c_pressure(external_node_id,:)=0;
    end
    Pi_g_p=PHI_NM_p*(PHI_MM_p\ones(size(PHI_MM_p,1),1));
    Pi_g_p=repmat(Pi_g_p,1,size(PHI_MM_p,1));
    Pi_g_p(Pi_g_p==0)=1;
    if interface_cube==0||interface_cube==2
        
        Pi_g_p_pressure=PHI_NM_p_pressure*(PHI_MM_p\ones(size(PHI_MM_p,1),1));
        Pi_g_p_pressure=repmat(Pi_g_p_pressure,1,size(PHI_MM_p,1));
        Pi_g_p_pressure(Pi_g_p_pressure==0)=1;
    end
    Pi_g_c=PHI_NM_c/PHI_MM_c*ones(size(PHI_MM_c,1),1);
    Pi_g_c=repmat(Pi_g_c,1,size(PHI_MM_c,1));
    Pi_g_c(Pi_g_c==0)=1;
    if interface_cube==0||interface_cube==2
        
        Pi_g_c_pressure=PHI_NM_c_pressure/PHI_MM_c*ones(size(PHI_MM_c,1),1);
        Pi_g_c_pressure=repmat(Pi_g_c_pressure,1,size(PHI_MM_c,1));
        Pi_g_c_pressure(Pi_g_c_pressure==0)=1;
    end
    Pr_node_cp=(PHI_NM_p/PHI_MM_p)./(Pi_g_p);
    Pr_node_pc=(PHI_NM_c/PHI_MM_c)./(Pi_g_c);
    if interface_cube==0||interface_cube==2
        
        Pr_node_cp_pressure=(PHI_NM_p_pressure/PHI_MM_p)./(Pi_g_p_pressure);
        Pr_node_pc_pressure=(PHI_NM_c_pressure/PHI_MM_c)./(Pi_g_c_pressure);
        Prpc_pressure=zeros(size(Mp,1),size(Mc,1));
    end
    Prpc=zeros(size(Mp,1),size(Mc,1));
    if interface_cube==0||interface_cube==2
        Prpc_pressure=zeros(size(Mc,1),size(Mp,1)).';
    end
    for nn=1:size(Pr_node_pc,1);
        for mm=1:size(Pr_node_pc,2)
            Prpc(3*((nn)-1)+(1:3),3*(mm-1)+(1:3)) = Pr_node_pc(nn,mm)*eye(3);
            if interface_cube==0||interface_cube==2
                Prpc_pressure(3*((nn)-1)+(1:3),3*(mm-1)+(1:3)) = Pr_node_pc_pressure(nn,mm)*eye(3);
            end
        end
    end
    Prcp=zeros(size(Mp,1),size(Mc,1)).';
    if interface_cube==0||interface_cube==2
        Prcp_pressure=zeros(size(Mp,1),size(Mc,1)).';
    end
    for nn=1:size(Pr_node_pc,1);
        for mm=1:size(Pr_node_pc,2)
            Prcp(3*(mm-1)+(1:3),3*((nn)-1)+(1:3)) = Pr_node_cp(mm,nn)*eye(3);
            if interface_cube==0||interface_cube==2
                Prcp_pressure(3*(mm-1)+(1:3),3*((nn)-1)+(1:3)) = Pr_node_cp_pressure(mm,nn)*eye(3);
            end
        end
    end
else
    ELintip=FACE_int_p;
    for k=1:size(FACE_int_p,1)
        for l=1:size(FACE_int_p,2)-1
            ELintip(k,l+1)=find(interface_p==FACE_int_p(k,l+1));
        end
    end
    ELintic=FACE_int_c;
    for k=1:size(FACE_int_c,1)
        for l=1:size(FACE_int_c,2)-1
            ELintic(k,l+1)=find(interface_c==FACE_int_c(k,l+1));
        end
    end
    [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,interface_p_coord);
    [csi_eta2]=identification_Q4function(Plate_coord,FACE_int_p,interface_c_coord);
    [Pr_node_pc]= Lagrange_based_interpolation(csi_eta,ELintic);
    [Pr_node_cp]= Lagrange_based_interpolation(csi_eta2,ELintip);
    %     Pr_node_pc(Pr_node_pc>1)=1;Pr_node_cp(Pr_node_cp>1)=1;
    %     Pr_node_pc(Pr_node_pc<0)=0;Pr_node_cp(Pr_node_cp<0)=0;
    if interface_cube==0||interface_cube==2
        Pr_node_pc_pressure=Pr_node_pc;
        Pr_node_pc_pressure(external_node_id,:)=0;
        Pr_node_cp_pressure=Pr_node_cp;
        Pr_node_cp_pressure(:,external_node_id)=0;
    end
    Prpc=zeros(size(Mp,1),size(Mc,1));
    if interface_cube==0||interface_cube==2
        Prpc_pressure=zeros(size(Mc,1),size(Mp,1)).';
    end
    for nn=1:size(Pr_node_pc,1);
        for mm=1:size(Pr_node_pc,2)
            Prpc(3*((nn)-1)+(1:3),3*(mm-1)+(1:3)) = Pr_node_pc(nn,mm)*eye(3);
            if interface_cube==0||interface_cube==2
                Prpc_pressure(3*((nn)-1)+(1:3),3*(mm-1)+(1:3)) = Pr_node_pc_pressure(nn,mm)*eye(3);
            end
        end
    end
    Prcp=zeros(size(Mp,1),size(Mc,1)).';
    if interface_cube==0||interface_cube==2
        Prcp_pressure=zeros(size(Mp,1),size(Mc,1)).';
    end
    for nn=1:size(Pr_node_pc,1);
        for mm=1:size(Pr_node_pc,2)
            Prcp(3*(mm-1)+(1:3),3*((nn)-1)+(1:3)) = Pr_node_cp(mm,nn)*eye(3);
            if interface_cube==0||interface_cube==2
                Prcp_pressure(3*(mm-1)+(1:3),3*((nn)-1)+(1:3)) = Pr_node_cp_pressure(mm,nn)*eye(3);
            end
        end
    end
end
%pressure filter
% When a node is outside the intersection of contact area the pressure is
% zero. The corrispondent line and column in Prcp and Prpc have to be
% put to zero.

[Ic,Jc,k_ijc,Isc,Jsc,DBsc]=Stiffness_assembly(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS);
if average_stress==1&&global_least_sq==0
    [Isc,Jsc,DBsc]=nodal_stress(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS);
elseif average_stress==1&&global_least_sq==1
    [Isc,Jsc,DBsc,Inc,Jnc,Nnnc]=nodal_stress(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS);
end
% [Ic,Jc,k_ijc,Isc,Jsc,DBsc]=nodal(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS);
Kc=sparse(Ic,Jc,k_ijc,3*length(Xc),3*length(Xc)); Kc=Kc+triu(Kc.',1);
[Ip,Jp,k_ijp,Isp,Jsp,DBsp]=Stiffness_assembly(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS);
if average_stress==1&&global_least_sq==0
    [Isp,Jsp,DBsp,]=nodal_stress(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS);
    elseif average_stress==1&&global_least_sq==1
        [Isp,Jsp,DBsp,Inp,Jnp,Nnnp]=nodal_stress(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS);
end
Kp=sparse(Ip,Jp,k_ijp,3*length(Xp),3*length(Xp)); Kp=Kp+triu(Kp.',1);
int_dof_plate=int_dofs_p(:,5);
int_dof_cube=int_dofs_c(:,5);
if Intergrid_switch==1
    Q_pc=(Mp*Prpc)/Mc; Q_cp=(Mc*Prcp)/Mp;
    if interface_cube==0||interface_cube==2
        Q_pc=(Mp*Prpc_pressure)/Mc; Q_cp=(Mc*Prcp_pressure)/Mp;
    end
    %     if interp_residual==1
    %         Q_pc=Mp*(Prcp)'/Mc;Q_cp=Mc*(Prpc)'/Mp;
    %     end
    if interp_residual==1
        
        Prcp=diag((Mc\ones(size(Mc,1),1))./(Prcp*(Mp\ones(size(Mp,1),1))))*Prcp;
        Prpc=diag((Mp\ones(size(Mp,1),1))./(Prpc*(Mc\ones(size(Mc,1),1))))*Prpc;
        Q_pc=Mp\(Prcp)'*Mc;Q_cp=Mc\(Prpc)'*Mp;
    end
    if correzione_momento==2
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        for i=1:size(Q_pc,2)
            A=[As;B(:,3*(i-1)+(1:3))'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            x0=Q_pc(:,i);
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Q_pc(:,i)=x;
        end
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        for i=1:size(Q_cp,2)
            A=[As2;B2(:,3*(i-1)+(1:3))'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            x0=Q_cp(:,i);
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Q_cp(:,i)=x;
        end
    end
elseif Intergrid_switch==0
    if apriori_patch==1
        r_sc=zeros(size(Mc,1),1);
        for k=1:size(Mc,1)/3
            r_sc(3*(k-1)+(1:3))=[0;0;sum(Mc(3*k,:))];
        end
        r_sp=zeros(size(Mp,1),1);
        for k=1:size(Mp,1)/3
            r_sp(3*(k-1)+(1:3))=[0;0;sum(Mp(3*k,:))];
        end
        r_sc1=zeros(size(Mc,1),1);
        for k=1:size(Mc,1)/3
            r_sc1(3*(k-1)+(1:3))=[sum(Mc(3*k-2,:));0;0];
        end
        r_sp1=zeros(size(Mp,1),1);
        for k=1:size(Mp,1)/3
            r_sp1(3*(k-1)+(1:3))=[sum(Mp(3*k-2,:));0;0];
        end
        r_sc2=zeros(size(Mc,1),1);
        for k=1:size(Mc,1)/3
            r_sc2(3*(k-1)+(1:3))=[0;sum(Mc(3*k-1,:));0];
        end
        r_sp2=zeros(size(Mp,1),1);
        for k=1:size(Mp,1)/3
            r_sp2(3*(k-1)+(1:3))=[0;sum(Mp(3*k,:));0];
        end
        Prcpn=Prcp*diag((sign(Prcp'*r_sc)+((Prcp'*r_sc)==0)).*max(abs((Prcp'*r_sc)),1e-6).\max(r_sp,1e-6));Prpcn=Prpc*diag((sign(Prpc'*r_sp)+((Prpc'*r_sp)==0)).*max(abs((Prpc'*r_sp)),1e-6).\max(r_sc,1e-6));
        Prcpn=Prcpn*diag((sign(Prcpn'*r_sc1)+((Prcpn'*r_sc1)==0)).*max(abs((Prcpn'*r_sc1)),1e-6).\max(r_sp1,1e-6));Prpcn=Prpcn*diag((sign(Prpcn'*r_sp1)+((Prpcn'*r_sp1)==0)).*max(abs((Prpcn'*r_sp1)),1e-6).\max(r_sc1,1e-6));
        Prcpn=Prcpn*diag((sign(Prcpn'*r_sc2)+((Prcpn'*r_sc2)==0)).*max(abs((Prcpn'*r_sc2)),1e-6).\max(r_sp2,1e-6));Prpcn=Prpcn*diag((sign(Prpcn'*r_sp2)+((Prpcn'*r_sp2)==0)).*max(abs((Prpcn'*r_sp2)),1e-6).\max(r_sc2,1e-6));
        if cube_master==0
            while max(max(abs(Prcp-Prcpn)))>1e-6
                Prcp=diag(sum(Prcpn,2))\Prcpn;
                Prcpn=Prcp*diag((sign(Prcp'*r_sc)+((Prcp'*r_sc)==0)).*max(abs((Prcp'*r_sc)),1e-6).\max(r_sp,1e-6));
                Prcpn=Prcpn*diag((sign(Prcpn'*r_sc1)+((Prcpn'*r_sc1)==0)).*max(abs((Prcpn'*r_sc1)),1e-6).\max(r_sp1,1e-6));
                Prcpn=Prcpn*diag((sign(Prcpn'*r_sc2)+((Prcpn'*r_sc2)==0)).*max(abs((Prcpn'*r_sc2)),1e-6).\max(r_sp2,1e-6));
            end
            
            Prcp=Prcpn;
        else
            iter=0;
            while max(max(abs(Prpc-Prpcn)))>1e-6
                iter=iter+1;
                Prpc=diag(sum(Prpcn,2))\Prpcn;
                Prpcn=Prpc*diag((sign(Prpc'*r_sp)+((Prpc'*r_sp)==0)).*max(abs((Prpc'*r_sp)),1e-6).\max(r_sc,1e-6));
                Prpcn=Prpcn*diag((sign(Prpcn'*r_sp1)+((Prpcn'*r_sp1)==0)).*max(abs((Prpcn'*r_sp1)),1e-6).\max(r_sc1,1e-6));
                Prcpn=Prcpn*diag((sign(Prcpn'*r_sc2)+((Prcpn'*r_sc2)==0)).*max(abs((Prcpn'*r_sc2)),1e-6).\max(r_sp2,1e-6));
                scatter(iter, max(max(abs(Prpc-Prpcn))),'fill')
                drawnow
                hold on
            end
            Prpc=Prpcn;
        end
    end
    Q_pc=Prcp';Q_cp=Prpc';
    if interp_residual==1
        Q_pc=Mp\(Prcp)'*Mc;Q_cp=Mc\(Prpc)'*Mp;
        Prcp=diag(Mc\ones(size(Mc,1),1)./(Prcp*(Mp\ones(size(Mp,1),1))))\Prcp;
        Prpc=diag(Mp\ones(size(Mp,1),1)./(Prpc*(Mc\ones(size(Mc,1),1))))\Prpc;
    end
elseif Intergrid_switch==2
    Q_pc=Prpc_pressure./repmat(sum(Prpc_pressure),size(Prpc_pressure,1),1);Q_cp=Prcp_pressure./repmat(sum(Prcp_pressure),size(Prcp_pressure,1),1);
    %     if interp_residual==1
    %         Q_pc=Mp*(Prcp)'/Mc;Q_cp=Mc*(Prpc)'/Mp;
    %     end
elseif Intergrid_switch==3
    Q_pc=Mp*Prpc_pressure/Mc./repmat(sum(Mp*Prpc_pressure/Mc),size(Mp*Prpc_pressure/Mc,1),1);Q_cp=Mc*Prcp_pressure/Mp./repmat(sum(Mc*Prcp_pressure/Mp),size(Mc*Prcp_pressure/Mp,1),1);
elseif Intergrid_switch==4
    Q_pc=Prpc_pressure./repmat(sum(Mp*Prpc_pressure/Mc),size(Mp*Prpc_pressure/Mc,1),1);Q_cp=Prcp_pressure./repmat(sum(Mc*Prcp_pressure/Mp),size(Mc*Prcp_pressure/Mp,1),1);
    Q_pc=Mp*Q_pc/Mc;Q_cp=Mc*Q_cp/Mp;
elseif Intergrid_switch==5
    Q_pc=(((eye(size(Prpc,2))+Prpc'*Prpc))\(Prcp+Prpc'))';Q_cp=((eye(size(Prcp,2))+Prcp'*Prcp)\(Prpc+Prcp'))';
    Prcp=((eye(size(Prpc,2))+Prpc'*Prpc))\(Prcp+Prpc');Prpc=((eye(size(Prcp,2))+Prcp'*Prcp)\(Prpc+Prcp'));
elseif Intergrid_switch==6
    [GP_coord,El,detJ]=Gauss_point(FACE_int_p,Plate_coord);
    [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,GP_coord);
    [Mcp]=intergrid_mass_matrix(csi_eta,El,detJ,FACE_int_c,FACE_int_p,Cube_DOFs_Map,Plate_DOFs_Map);
    Mcp=Mcp(int_dofs_c(:,5),int_dofs_p(:,5));
    [GP_coord,El,detJ]=Gauss_point(FACE_int_c,Cube_coord);
    [csi_eta]=identification_Q4function(Plate_coord,FACE_int_p,GP_coord);
    [Mpc]=intergrid_mass_matrix(csi_eta,El,detJ,FACE_int_p,FACE_int_c,Plate_DOFs_Map,Cube_DOFs_Map);
    Mpc=Mpc(int_dofs_p(:,5),int_dofs_c(:,5));
    Q_pc=Mpc/Mc; Q_cp=Mcp/Mp;
elseif Intergrid_switch==7
    %Mortar Symmetric
    [GP_coord,El,detJ]=Gauss_point(FACE_int_p,Plate_coord);
    [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,GP_coord);
    [Mcp]=intergrid_mass_matrix(csi_eta,El,detJ,FACE_int_c,FACE_int_p,Cube_DOFs_Map,Plate_DOFs_Map);
    Mcp=Mcp(int_dofs_c(:,5),int_dofs_p(:,5));
    [GP_coord,El,detJ]=Gauss_point(FACE_int_c,Cube_coord);
    [csi_eta2]=identification_Q4function(Plate_coord,FACE_int_p,GP_coord);
    [Mpc]=intergrid_mass_matrix(csi_eta2,El,detJ,FACE_int_p,FACE_int_c,Plate_DOFs_Map,Cube_DOFs_Map);
    Mpc=Mpc(int_dofs_p(:,5),int_dofs_c(:,5));
    Prpc=Mp\Mpc;Prcp=Mc\Mcp;
    Q_pc=Prcp.';Q_cp=Prpc.';
elseif Intergrid_switch==8
    %Internodes corrigé
    %     [GP_coord,El,detJ]=Gauss_point(FACE_int_p,Plate_coord);
    [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,interface_p_coord);
    [csi_eta2]=identification_Q4function(Plate_coord,FACE_int_p,interface_c_coord);
    [Mortar_GPpc,mortar_detjpc]=Geometric_intersection_quadratic(interface_c_coord,ELintic,interface_p_coord,ELintip,csi_eta,csi_eta2);
    [Mortar_GPcp,mortar_detjcp]=Geometric_intersection_quadratic(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    [Mc]=intergrid_mass_matrix_intersection(Mortar_GPpc,mortar_detjpc,FACE_int_c,FACE_int_p,Cube_DOFs_Map);
    
    Mc=Mc(int_dofs_c(:,5),int_dofs_c(:,5));
    
    %     [Mortar_GP,mortar_detj]=Geometric_intersection_quadratic(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    [Mp]=intergrid_mass_matrix_intersection(Mortar_GPcp,mortar_detjcp,FACE_int_p,FACE_int_c,Plate_DOFs_Map);
    Mp=Mp(int_dofs_p(:,5),int_dofs_p(:,5));
    if interface_cube==0||interface_cube==2
        Q_pc=Mp*Prpc_pressure/Mc; Q_cp=Mc*Prcp_pressure/Mp;
    else
        Q_pc=Mp*Prpc/Mc; Q_cp=Mc*Prcp/Mp;
    end
elseif Intergrid_switch==9
    %Internodes corrigé et rescalé
    %     [GP_coord,El,detJ]=Gauss_point(FACE_int_p,Plate_coord);
    [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,interface_p_coord);
    [csi_eta2]=identification_Q4function(Plate_coord,FACE_int_p,interface_c_coord);
    [Mortar_GPpc,mortar_detjpc]=Geometric_intersection_quadratic(interface_c_coord,ELintic,interface_p_coord,ELintip,csi_eta,csi_eta2);
    [Mortar_GPcp,mortar_detjcp]=Geometric_intersection_quadratic(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    [Mc]=intergrid_mass_matrix_intersection(Mortar_GP,mortar_detj,FACE_int_c,FACE_int_p,Cube_DOFs_Map);
    Mc=Mc(int_dofs_c(:,5),int_dofs_c(:,5));
    
    %     [Mortar_GP,mortar_detj]=Geometric_intersection(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    [Mp]=intergrid_mass_matrix_intersection(Mortar_GP,mortar_detj,FACE_int_p,FACE_int_c,Plate_DOFs_Map);
    Mp=Mp(int_dofs_p(:,5),int_dofs_p(:,5));
    
    if interface_cube==0||interface_cube==2
        Q_pc=Mp*Prpc_pressure/Mc./repmat(sum(Mp*Prpc_pressure/Mc),size(Mp*Prpc_pressure/Mc,1),1);Q_cp=Mc*Prcp_pressure/Mp./repmat(sum(Mc*Prcp_pressure/Mp),size(Mc*Prcp_pressure/Mp,1),1);
    else
        Q_pc=Mp*Prpc/Mc./repmat(sum(Mp*Prpc/Mc),size(Mp*Prpc/Mc,1),1);Q_cp=Mc*Prcp/Mp./repmat(sum(Mc*Prcp/Mp),size(Mc*Prcp/Mp,1),1);
    end
    
elseif Intergrid_switch==10
    %Mortar corrigé
    %     [GP_coord,El,detJ]=Gauss_point(FACE_int_p,Plate_coord);
    [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,interface_p_coord);
    [csi_eta2]=identification_Q4function(Plate_coord,FACE_int_p,interface_c_coord);
    [Mortar_GPpc,mortar_detjpc]=Geometric_intersection_quadratic(interface_c_coord,ELintic,interface_p_coord,ELintip,csi_eta,csi_eta2);
    [Mortar_GPcp,mortar_detjcp]=Geometric_intersection_quadratic(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    %     [Mc]=intergrid_mass_matrix_intersection(Mortar_GPcp,mortar_detjcp,FACE_int_c,FACE_int_p,Cube_DOFs_Map);
    [Mc]=intergrid_mass_matrix_intersection(Mortar_GPpc,mortar_detjpc,FACE_int_c,FACE_int_p,Cube_DOFs_Map);
    
    Mc=Mc(int_dofs_c(:,5),int_dofs_c(:,5));
    %     [Mp]=intergrid_mass_matrix_intersection(Mortar_GPpc,mortar_detjpc,FACE_int_p,FACE_int_c,Plate_DOFs_Map);
    [Mp]=intergrid_mass_matrix_intersection(Mortar_GPcp,mortar_detjcp,FACE_int_p,FACE_int_c,Plate_DOFs_Map);
    
    Mp=Mp(int_dofs_p(:,5),int_dofs_p(:,5));
    [Mpc]=intergrid_mass_matrix_special(Mortar_GPcp,Mortar_GPpc,mortar_detjcp,FACE_int_p,FACE_int_c,Cube_DOFs_Map,Plate_DOFs_Map);
    [Mcp]=intergrid_mass_matrix_special(Mortar_GPpc,Mortar_GPcp,mortar_detjpc,FACE_int_c,FACE_int_p,Plate_DOFs_Map,Cube_DOFs_Map);
    Mcp=Mcp(int_dofs_c(:,5),int_dofs_p(:,5));
    Mpc=Mpc(int_dofs_p(:,5),int_dofs_c(:,5));
    Prpc=Mp\Mpc;Prcp=Mc\Mcp;
    Q_pc=Prcp.';Q_cp=Prpc.';
    if correzione_momento==1
        %tensore dei bracci.
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        for i=1:size(Prcp,1)
            A=repmat((Prcp(i,:)~=0)',1,3);
            C=(null(full(B(:,3*(i-1)+(1:3)).*A)'))';
            Prcp(i,:)=(Prcp(i,:)*C.')*C;
            Prcp(i,:)=Prcp(i,:)/(sum(Prcp(i,:))*(sum(Prcp(i,:))~=0)+1*(sum(Prcp(i,:))==0));
        end
        for i=1:size(Prpc,1)
            A=repmat((Prpc(i,:)~=0)',1,3);
            C=(null(full(B2(:,3*(i-1)+(1:3)).*A)'))';
            Prpc(i,:)=(Prpc(i,:)*C.')*C;
            Prpc(i,:)=Prpc(i,:)/(sum(Prpc(i,:))*(sum(Prpc(i,:))~=0)+1*(sum(Prpc(i,:))==0));
        end
        %         C=(null((B)'))';
        %         C2=(null((B2)'))';
        %         Prcp=(Prcp*C.')*C;Prpc=(Prpc*C2.')*C2;
        %         Prcp=diag(sum(Prcp,2))\Prcp;
        %         Prpc=diag(sum(Prpc,2))\Prpc;
        
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==2
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        for i=1:size(Prcp,1)
            A=[As;B(:,3*(i-1)+(1:3))'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            x0=Prcp(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prcp(i,:)=x';
        end
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        for i=1:size(Prpc,1)
            A=[As2;B2(:,3*(i-1)+(1:3))'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            x0=Prpc(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prpc(i,:)=x';
        end
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==3
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        %analytic solution constant stress
        for nc=1:size(interface_c_coord,1)
            Dc(3*(nc-1)+(1:3),1:3)=repmat(interface_c_coord(nc,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1]+1;
        end
        for np=1:size(interface_p_coord,1)
            Dp(3*(np-1)+(1:3),1:3)=repmat(interface_p_coord(np,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1]+1;
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        for i=1:size(Prcp,1)
            A=[As;B(:,3*(i-1)+(1:3))';Dp'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            b=[b;Dc(i,:)'];
            x0=Prcp(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prcp(i,:)=x';
        end
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        for i=1:size(Prpc,1)
            A=[As2;B2(:,3*(i-1)+(1:3))';Dc'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            b=[b;Dp(i,:)'];
            x0=Prpc(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prpc(i,:)=x';
        end
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==4
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        %analytic solution constant stress
        for nc=1:size(interface_c_coord,1)
            Dc(3*(nc-1)+(1:3),1:3)=repmat(interface_c_coord(nc,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1];
        end
        for np=1:size(interface_p_coord,1)
            Dp(3*(np-1)+(1:3),1:3)=repmat(interface_p_coord(np,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1];
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        Rsc=As2*Mc; Rsp=As*Mp;
        [Icp,Jcp]=find(ones(size(Prcp)));
        [Ipc,Jpc]=find(ones(size(Prpc)));
        A=zeros(9*size(Prcp,1)+3*size(Prcp,2),length(Prcp(:)));
        Bs=zeros(9*size(Prcp,1)+3*size(Prcp,2),1);
        for i=1:size(Prcp,1)
            A((i-1)*9+(1:9),Icp==i)=[As;B(:,3*(i-1)+(1:3))';Dp'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            Bs((i-1)*9+(1:9))=[b;Dc(i,:)'];
        end
        for j=1:size(Prcp,2)
            A(9*size(Prcp,1)+(j-1)*3+(1:3),Jcp==j)=Rsc;
            Bs(9*size(Prcp,1)+(j-1)*3+(1:3))=Rsp(:,j);
        end
        A=sparse(A);
        Bs=sparse(Bs);
        x0=Prcp(:);
        x0=sparse(x0);
        Lambda=2*((A*A')\(Bs-A*x0));
        x=x0+A'*Lambda/2;
        Prcp=reshape(x,size(Prcp,1),size(Prcp,2));
        A=zeros(9*size(Prcp,1)+3*size(Prcp,2),length(Prcp(:)));
        Bs=zeros(9*size(Prcp,1)+3*size(Prcp,2),1);
        for i=1:size(Prpc,1)
            A((i-1)*9+(1:9),Icp==i)=[As2;B2(:,3*(i-1)+(1:3))';Dc'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            Bs((i-1)*9+(1:9))=[b;Dp(i,:)'];
        end
        for j=1:size(Prpc,2)
            A(9*size(Prpc,1)+(j-1)*3+(1:3),Jpc==j)=Rsp;
            Bs(9*size(Prpc,1)+(j-1)*3+(1:3))=Rsc(:,j);
        end
        A=sparse(A);
        Bs=sparse(Bs);
        x0=Prpc(:);
        x0=sparse(x0);
        Lambda=2*((A*A')\(Bs-A*x0));
        x=x0+A'*Lambda/2;
        Prpc=reshape(x,size(Prpc,1),size(Prpc,2));
        Q_cp=Prpc.'; Q_pc=Prcp.';
    end
    if interp_residual==1
        Q_pc=Mp\(Prcp)'*Mc;Q_cp=Mc\(Prpc)'*Mp;
    end
elseif Intergrid_switch==11
    %Bilan local
    %     [GP_coord,El,detJ]=Gauss_point(FACE_int_p,Plate_coord);
    [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,interface_p_coord);
    [csi_eta2]=identification_Q4function(Plate_coord,FACE_int_p,interface_c_coord);
    [Mortar_GPpc,mortar_detjpc]=Geometric_intersection_quadratic(interface_c_coord,ELintic,interface_p_coord,ELintip,csi_eta,csi_eta2);
    [Mortar_GPcp,mortar_detjcp]=Geometric_intersection_quadratic(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    [Mc]=intergrid_mass_matrix_intersection(Mortar_GPpc,mortar_detjpc,FACE_int_c,FACE_int_p,Cube_DOFs_Map);
    Mc=Mc(int_dofs_c(:,5),int_dofs_c(:,5));
    [Mp]=intergrid_mass_matrix_intersection(Mortar_GPcp,mortar_detjcp,FACE_int_p,FACE_int_c,Plate_DOFs_Map);
    Mp=Mp(int_dofs_p(:,5),int_dofs_p(:,5));
    [Mpc]=intergrid_mass_matrix_special(Mortar_GPcp,Mortar_GPpc,mortar_detjpc,FACE_int_p,FACE_int_c,Cube_DOFs_Map,Plate_DOFs_Map);
    [Mcp]=intergrid_mass_matrix_special(Mortar_GPpc,Mortar_GPcp,mortar_detjcp,FACE_int_c,FACE_int_p,Plate_DOFs_Map,Cube_DOFs_Map);
    Mcp=Mcp(int_dofs_c(:,5),int_dofs_p(:,5));
    Mpc=Mpc(int_dofs_p(:,5),int_dofs_c(:,5));
    Q_pc=Mpc/Mc; Q_cp=Mcp/Mp;
elseif Intergrid_switch==12
    %Simone 1 (not enough)
    %     [GP_coord,El,detJ]=Gauss_point(FACE_int_p,Plate_coord);
    [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,interface_p_coord);
    [csi_eta2]=identification_Q4function(Plate_coord,FACE_int_p,interface_c_coord);
    [Mortar_GPpc,mortar_detjpc]=Geometric_intersection_quadratic(interface_c_coord,ELintic,interface_p_coord,ELintip,csi_eta,csi_eta2);
    [Mortar_GPcp,mortar_detjcp]=Geometric_intersection_quadratic(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    [kpc,ikp,jkc]=normal_stress_faces(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS,FACE_int_c,Mortar_GPcp,mortar_detjcp,Mortar_GPpc,FACE_int_p);
    Kpc=sparse(ikp,jkc,kpc,max(max(Plate_DOFs_Map)),max(max(CUBE_ZONE_ELEMENT_DOFS)));
    Kpc=Kpc(int_dofs_p(:,5),:);
    [kcp,ikc,jkp]=normal_stress_faces(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS,FACE_int_p,Mortar_GPpc,mortar_detjpc,Mortar_GPcp,FACE_int_c);
    Kcp=sparse(ikc,jkp,kcp,max(max(Cube_DOFs_Map)),max(max(PLATE_ZONE_ELEMENT_DOFS)));
    Kcp=Kcp(int_dofs_c(:,5),:);
elseif Intergrid_switch==13
    %Simone 2 (Lagrange eigenvalue based)
    %     [GP_coord,El,detJ]=Gauss_point(FACE_int_p,Plate_coord);
    [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,interface_p_coord);
    [csi_eta2]=identification_Q4function(Plate_coord,FACE_int_p,interface_c_coord);
    [Mortar_GPpc,mortar_detjpc]=Geometric_intersection_quadratic(interface_c_coord,ELintic,interface_p_coord,ELintip,csi_eta,csi_eta2);
    [Mortar_GPcp,mortar_detjcp]=Geometric_intersection_quadratic(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    [kpc,ikp,jkc]=normal_stress_faces(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS,FACE_int_c,Mortar_GPcp,mortar_detjpc,Mortar_GPpc,FACE_int_p);
    Kpc=sparse(ikp,jkc,kpc,max(max(Plate_DOFs_Map)),max(max(CUBE_ZONE_ELEMENT_DOFS)));
    Kpc=Kpc(int_dofs_p(:,5),:);
    [kcp,ikc,jkp]=normal_stress_faces(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS,FACE_int_p,Mortar_GPpc,mortar_detjcp,Mortar_GPcp,FACE_int_c);
    Kcp=sparse(ikc,jkp,kcp,max(max(Cube_DOFs_Map)),max(max(PLATE_ZONE_ELEMENT_DOFS)));
    Kcp=Kcp(int_dofs_c(:,5),:);
    [Kcic]=normal_stress_faces_auto(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS,FACE_int_c,mortar_detjpc,Mortar_GPpc);
    Kcic=Kcic(int_dofs_c(:,5),:);
    [Kpip]=normal_stress_faces_auto(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS,FACE_int_p,mortar_detjcp,Mortar_GPcp);
    Kpip=Kpip(int_dofs_p(:,5),:);
    %add the displacement decomposition to avoid rigid mode
    %evaluate the rigid mode matrix of the interface
    %3 translations per node and
    Rigid_c=repmat(speye(3),size(interface_c_coord,1),1);
    Rigid_c=[Rigid_c,reshape(cross(repmat([1;0;0],1,size(interface_c_coord,1)),interface_c_coord'),[],1),reshape(cross(repmat([0;1;0],1,size(interface_c_coord,1)),interface_c_coord'),[],1),reshape(cross(repmat([0;0;1],1,size(interface_c_coord,1)),interface_c_coord'),[],1)];
    Rigid_p=repmat(speye(3),size(interface_p_coord,1),1);
    Rigid_p=[Rigid_p,reshape(cross(repmat([1;0;0],1,size(interface_p_coord,1)),interface_p_coord'),[],1),reshape(cross(repmat([0;1;0],1,size(interface_p_coord,1)),interface_p_coord'),[],1),reshape(cross(repmat([0;0;1],1,size(interface_p_coord,1)),interface_p_coord'),[],1)];
    Defc=null(full(Rigid_c)');
    Defp=null(full(Rigid_p)');
elseif Intergrid_switch==14
    %Simone 3 (Lagrange eigenvalue based with stress flux and displacement jump)
    %     [GP_coord,El,detJ]=Gauss_point(FACE_int_p,Plate_coord);
    [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,interface_p_coord);
    [csi_eta2]=identification_Q4function(Plate_coord,FACE_int_p,interface_c_coord);
    [Mortar_GPpc,mortar_detjpc]=Geometric_intersection_quadratic(interface_c_coord,ELintic,interface_p_coord,ELintip,csi_eta,csi_eta2);
    [Mortar_GPcp,mortar_detjcp]=Geometric_intersection_quadratic(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    [kpc,ikp,jkc]=normal_stress_faces(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS,FACE_int_c,Mortar_GPcp,mortar_detjpc,Mortar_GPpc,FACE_int_p);
    Kpc=sparse(ikp,jkc,kpc,max(max(Plate_DOFs_Map)),max(max(CUBE_ZONE_ELEMENT_DOFS)));
    Kpc=Kpc(int_dofs_p(:,5),:);
    [kcp,ikc,jkp]=normal_stress_faces(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS,FACE_int_p,Mortar_GPpc,mortar_detjcp,Mortar_GPcp,FACE_int_c);
    Kcp=sparse(ikc,jkp,kcp,max(max(Cube_DOFs_Map)),max(max(PLATE_ZONE_ELEMENT_DOFS)));
    Kcp=Kcp(int_dofs_c(:,5),:);
    [Kcic]=normal_stress_faces_auto(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS,FACE_int_c,mortar_detjpc,Mortar_GPpc);
    Kcic=Kcic(int_dofs_c(:,5),:);
    [Kpip]=normal_stress_faces_auto(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS,FACE_int_p,mortar_detjcp,Mortar_GPcp);
    Kpip=Kpip(int_dofs_p(:,5),:);
    [Mc]=intergrid_mass_matrix_intersection(Mortar_GPpc,mortar_detjpc,FACE_int_c,FACE_int_p,Cube_DOFs_Map);
    
    Mc=Mc(int_dofs_c(:,5),int_dofs_c(:,5));
    %     [Mp]=intergrid_mass_matrix_intersection(Mortar_GPpc,mortar_detjpc,FACE_int_p,FACE_int_c,Plate_DOFs_Map);
    [Mp]=intergrid_mass_matrix_intersection(Mortar_GPcp,mortar_detjcp,FACE_int_p,FACE_int_c,Plate_DOFs_Map);
    
    Mp=Mp(int_dofs_p(:,5),int_dofs_p(:,5));
    [Mpc]=intergrid_mass_matrix_special(Mortar_GPcp,Mortar_GPpc,mortar_detjcp,FACE_int_p,FACE_int_c,Cube_DOFs_Map,Plate_DOFs_Map);
    [Mcp]=intergrid_mass_matrix_special(Mortar_GPpc,Mortar_GPcp,mortar_detjpc,FACE_int_c,FACE_int_p,Plate_DOFs_Map,Cube_DOFs_Map);
    Mcp=Mcp(int_dofs_c(:,5),int_dofs_p(:,5));
    Mpc=Mpc(int_dofs_p(:,5),int_dofs_c(:,5));
    %add the displacement decomposition to avoid rigid mode
    %evaluate the rigid mode matrix of the interface
    %3 translations per node and
    %     Rigid_c=repmat(speye(3),size(interface_c_coord,1),1);
    %     Rigid_c=[Rigid_c,reshape(cross(repmat([1;0;0],1,size(interface_c_coord,1)),interface_c_coord'),[],1),reshape(cross(repmat([0;1;0],1,size(interface_c_coord,1)),interface_c_coord'),[],1),reshape(cross(repmat([0;0;1],1,size(interface_c_coord,1)),interface_c_coord'),[],1)];
    %     Rigid_p=repmat(speye(3),size(interface_p_coord,1),1);
    %     Rigid_p=[Rigid_p,reshape(cross(repmat([1;0;0],1,size(interface_p_coord,1)),interface_p_coord'),[],1),reshape(cross(repmat([0;1;0],1,size(interface_p_coord,1)),interface_p_coord'),[],1),reshape(cross(repmat([0;0;1],1,size(interface_p_coord,1)),interface_p_coord'),[],1)];
    %     Defc=null(full(Rigid_c)');
    %     Defp=null(full(Rigid_p)');
elseif Intergrid_switch==15
    %Simone 4 (Lagrange eigenvalue based with all stress flux and displacement jump)
    %     [GP_coord,El,detJ]=Gauss_point(FACE_int_p,Plate_coord);
    [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,interface_p_coord);
    [csi_eta2]=identification_Q4function(Plate_coord,FACE_int_p,interface_c_coord);
    [Mortar_GPpc,mortar_detjpc]=Geometric_intersection_quadratic(interface_c_coord,ELintic,interface_p_coord,ELintip,csi_eta,csi_eta2);
    [Mortar_GPcp,mortar_detjcp]=Geometric_intersection_quadratic(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    [kpc,ikp,jkc]=normal_stress_faces(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS,FACE_int_c,Mortar_GPcp,mortar_detjpc,Mortar_GPpc,FACE_int_p);
    Kpc=sparse(ikp,jkc,kpc,max(max(Plate_DOFs_Map)),max(max(CUBE_ZONE_ELEMENT_DOFS)));
    Kpc=Kpc(int_dofs_p(:,5),:);
    [kcp,ikc,jkp]=normal_stress_faces(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS,FACE_int_p,Mortar_GPpc,mortar_detjcp,Mortar_GPcp,FACE_int_c);
    Kcp=sparse(ikc,jkp,kcp,max(max(Cube_DOFs_Map)),max(max(PLATE_ZONE_ELEMENT_DOFS)));
    Kcp=Kcp(int_dofs_c(:,5),:);
    [Kcic]=normal_stress_faces_auto(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS,FACE_int_c,mortar_detjpc,Mortar_GPpc);
    Kcic=Kcic(int_dofs_c(:,5),:);
    [Kpip]=normal_stress_faces_auto(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS,FACE_int_p,mortar_detjcp,Mortar_GPcp);
    Kpip=Kpip(int_dofs_p(:,5),:);
    %tangential
    [kpc,ikp,jkc]=tangential_stress_faces(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS,FACE_int_c,Mortar_GPcp,mortar_detjpc,Mortar_GPpc,FACE_int_p);
    Ktpc=sparse(ikp,jkc,kpc,max(max(Plate_DOFs_Map)),max(max(CUBE_ZONE_ELEMENT_DOFS)));
    Ktpc=Ktpc(int_dofs_p(:,5),:);
    [kcp,ikc,jkp]=tangential_stress_faces(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS,FACE_int_p,Mortar_GPpc,mortar_detjcp,Mortar_GPcp,FACE_int_c);
    Ktcp=sparse(ikc,jkp,kcp,max(max(Cube_DOFs_Map)),max(max(PLATE_ZONE_ELEMENT_DOFS)));
    Ktcp=Ktcp(int_dofs_c(:,5),:);
    [Ktcic]=tangential_stress_faces_auto(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS,FACE_int_c,mortar_detjpc,Mortar_GPpc);
    Ktcic=Ktcic(int_dofs_c(:,5),:);
    [Ktpip]=tangential_stress_faces_auto(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS,FACE_int_p,mortar_detjcp,Mortar_GPcp);
    Ktpip=Ktpip(int_dofs_p(:,5),:);
    %bitan
    [kpc,ikp,jkc]=bitangential_stress_faces(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS,FACE_int_c,Mortar_GPcp,mortar_detjpc,Mortar_GPpc,FACE_int_p);
    Kbpc=sparse(ikp,jkc,kpc,max(max(Plate_DOFs_Map)),max(max(CUBE_ZONE_ELEMENT_DOFS)));
    Kbpc=Kbpc(int_dofs_p(:,5),:);
    [kcp,ikc,jkp]=bitangential_stress_faces(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS,FACE_int_p,Mortar_GPpc,mortar_detjcp,Mortar_GPcp,FACE_int_c);
    Kbcp=sparse(ikc,jkp,kcp,max(max(Cube_DOFs_Map)),max(max(PLATE_ZONE_ELEMENT_DOFS)));
    Kbcp=Kbcp(int_dofs_c(:,5),:);
    [Kbcic]=bitangential_stress_faces_auto(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS,FACE_int_c,mortar_detjpc,Mortar_GPpc);
    Kbcic=Kbcic(int_dofs_c(:,5),:);
    [Kbpip]=bitangential_stress_faces_auto(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS,FACE_int_p,mortar_detjcp,Mortar_GPcp);
    Kbpip=Kbpip(int_dofs_p(:,5),:);
    %mass
    [Mc]=intergrid_mass_matrix_intersection(Mortar_GPpc,mortar_detjpc,FACE_int_c,FACE_int_p,Cube_DOFs_Map);
    
    Mc=Mc(int_dofs_c(:,5),int_dofs_c(:,5));
    %     [Mp]=intergrid_mass_matrix_intersection(Mortar_GPpc,mortar_detjpc,FACE_int_p,FACE_int_c,Plate_DOFs_Map);
    [Mp]=intergrid_mass_matrix_intersection(Mortar_GPcp,mortar_detjcp,FACE_int_p,FACE_int_c,Plate_DOFs_Map);
    
    Mp=Mp(int_dofs_p(:,5),int_dofs_p(:,5));
    [Mpc]=intergrid_mass_matrix_special(Mortar_GPcp,Mortar_GPpc,mortar_detjcp,FACE_int_p,FACE_int_c,Cube_DOFs_Map,Plate_DOFs_Map);
    [Mcp]=intergrid_mass_matrix_special(Mortar_GPpc,Mortar_GPcp,mortar_detjpc,FACE_int_c,FACE_int_p,Plate_DOFs_Map,Cube_DOFs_Map);
    Mcp=Mcp(int_dofs_c(:,5),int_dofs_p(:,5));
    Mpc=Mpc(int_dofs_p(:,5),int_dofs_c(:,5));
    %add the displacement decomposition to avoid rigid mode
    %evaluate the rigid mode matrix of the interface
    %3 translations per node and
    %     Rigid_c=repmat(speye(3),size(interface_c_coord,1),1);
    %     Rigid_c=[Rigid_c,reshape(cross(repmat([1;0;0],1,size(interface_c_coord,1)),interface_c_coord'),[],1),reshape(cross(repmat([0;1;0],1,size(interface_c_coord,1)),interface_c_coord'),[],1),reshape(cross(repmat([0;0;1],1,size(interface_c_coord,1)),interface_c_coord'),[],1)];
    %     Rigid_p=repmat(speye(3),size(interface_p_coord,1),1);
    %     Rigid_p=[Rigid_p,reshape(cross(repmat([1;0;0],1,size(interface_p_coord,1)),interface_p_coord'),[],1),reshape(cross(repmat([0;1;0],1,size(interface_p_coord,1)),interface_p_coord'),[],1),reshape(cross(repmat([0;0;1],1,size(interface_p_coord,1)),interface_p_coord'),[],1)];
    %     Defc=null(full(Rigid_c)');
    %     Defp=null(full(Rigid_p)');
elseif Intergrid_switch==16
    %Mortar approximé
    %     [GP_coord,El,detJ]=Gauss_point(FACE_int_p,Plate_coord);
    %     [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,interface_p_coord);
    %     [csi_eta2]=identification_Q4function(Plate_coord,FACE_int_p,interface_c_coord);
    %     [Mortar_GPpc,mortar_detjpc]=Geometric_intersection_quadratic(interface_c_coord,ELintic,interface_p_coord,ELintip,csi_eta,csi_eta2);
    %     [Mortar_GPcp,mortar_detjcp]=Geometric_intersection_quadratic(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    %     [Mc]=intergrid_mass_matrix_intersection(Mortar_GPpc,mortar_detjpc,FACE_int_c,FACE_int_p,Cube_DOFs_Map);
    %
    %     Mc=Mc(int_dofs_c(:,5),int_dofs_c(:,5));
    %
    %     %     [Mortar_GP,mortar_detj]=Geometric_intersection_quadratic(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    %     [Mp]=intergrid_mass_matrix_intersection(Mortar_GPcp,mortar_detjcp,FACE_int_p,FACE_int_c,Plate_DOFs_Map);
    %     Mp=Mp(int_dofs_p(:,5),int_dofs_p(:,5));
    if interface_cube==0||interface_cube==2
        %         Q_pc=Mp*Prpc_pressure/Mc; Q_cp=Mc*Prcp_pressure/Mp;
        Prpc=Mp\diag(sum(Mp,2))*Prpc*(diag(sum(Mc,2))\Mc); Prcp=Mc\diag(sum(Mc,2))*Prcp*(diag(sum(Mp,2))\Mp);
        Q_cp=Prpc.'; Q_pc=Prcp.';
        
    else
        Prpc=Mp\diag(sum(Mp,2))*Prpc*(diag(sum(Mc,2))\Mc); Prcp=Mc\diag(sum(Mc,2))*Prcp*(diag(sum(Mp,2))\Mp);
        Q_cp=Prpc.'; Q_pc=Prcp.';
        %         Q_pc=Mp*Prpc/Mc; Q_cp=Mc*Prcp/Mp;
    end
    if correzione_momento==1
        %tensore dei bracci.
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        for i=1:size(Prcp,1)
            A=repmat((Prcp(i,:)~=0)',1,3);
            C=(null((B(:,3*(i-1)+(1:3)).*A)'))';
            Prcp(i,:)=(Prcp(i,:)*C.')*C;
            Prcp(i,:)=Prcp(i,:)/(sum(Prcp(i,:))*(sum(Prcp(i,:))~=0)+1*(sum(Prcp(i,:))==0));
        end
        for i=1:size(Prpc,1)
            A=repmat((Prpc(i,:)~=0)',1,3);
            C=(null((B2(:,3*(i-1)+(1:3)).*A)'))';
            Prpc(i,:)=(Prpc(i,:)*C.')*C;
            Prpc(i,:)=Prpc(i,:)/(sum(Prpc(i,:))*(sum(Prpc(i,:))~=0)+1*(sum(Prpc(i,:))==0));
        end
        %         C=(null((B)'))';
        %         C2=(null((B2)'))';
        %         Prcp=(Prcp*C.')*C;Prpc=(Prpc*C2.')*C2;
        %         Prcp=diag(sum(Prcp,2))\Prcp;
        %         Prpc=diag(sum(Prpc,2))\Prpc;
        
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==2
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        for i=1:size(Prcp,1)
            A=[As;B(:,3*(i-1)+(1:3))'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            x0=Prcp(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prcp(i,:)=x';
        end
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        for i=1:size(Prpc,1)
            A=[As2;B2(:,3*(i-1)+(1:3))'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            x0=Prpc(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prpc(i,:)=x';
        end
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==3
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        %analytic solution constant stress
        for nc=1:size(interface_c_coord,1)
            Dc(3*(nc-1)+(1:3),1)=repmat(interface_c_coord(nc,:),1,1).*[-.29 -.29 1];
        end
        for np=1:size(interface_p_coord,1)
            Dp(3*(np-1)+(1:3),1)=repmat(interface_p_coord(np,:),1,1).*[-.29 -.29 1];
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        for i=1:size(Prcp,1)
            A=[As;B(:,3*(i-1)+(1:3))';Dp'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            b=[b;Dc(i,:)'];
            x0=Prcp(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prcp(i,:)=x';
        end
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        for i=1:size(Prpc,1)
            A=[As2;B2(:,3*(i-1)+(1:3))';Dc'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            b=[b;Dp(i,:)'];
            x0=Prpc(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prpc(i,:)=x';
        end
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==4
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        %analytic solution constant stress
        %         for nc=1:size(interface_c_coord,1)
        %             Dc(3*(nc-1)+(1:3),1:3)=repmat(interface_c_coord(nc,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1]+10;
        %         end
        %         for np=1:size(interface_p_coord,1)
        %             Dp(3*(np-1)+(1:3),1:3)=repmat(interface_p_coord(np,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1]+10;
        %         end
        for nc=1:size(interface_c_coord,1)
            Dc(3*(nc-1)+(1:3),1)=repmat(interface_c_coord(nc,:),1,1).*[-.29 -.29 1];
        end
        for np=1:size(interface_p_coord,1)
            Dp(3*(np-1)+(1:3),1)=repmat(interface_p_coord(np,:),1,1).*[-.29 -.29 1];
        end
        for nc=1:size(Cube_coord,1)
            Ucs(3*(nc-1)+(1:3),1)=repmat(Cube_coord(nc,:),1,1).*[-.29 -.29 1];
        end
        for np=1:size(Plate_coord,1)
            Ups(3*(np-1)+(1:3),1)=repmat(Plate_coord(np,:),1,1).*[-.29 -.29 1];
        end
        Rcs=Kc*Ucs;Rps=Kp*Ups;
        As=repmat(eye(3),1,size(Prcp,2)/3);
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        %         Rsc=As2(3,:)*Mc; Rsp=As(3,:)*Mp;
        Rsc=Rcs(int_dof_cube)'; Rsp=Rps(int_dof_plate)';
        [Icp,Jcp]=find(ones(size(Prcp)));
        [Ipc,Jpc]=find(ones(size(Prpc)));
        A=zeros(7*size(Prcp,1)+size(Prcp,2),length(Prcp(:)));
        Bs=zeros(7*size(Prcp,1)+size(Prcp,2),1);
        %         A=zeros(6*size(Prcp,1),length(Prcp(:)));
        %         Bs=zeros(6*size(Prcp,1),1);
        for i=1:size(Prcp,1)
            A((i-1)*7+(1:7),Icp==i)=[As;B(:,3*(i-1)+(1:3))';Dp'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            Bs((i-1)*7+(1:7))=[b;Dc(i,:)'];
        end
        for j=1:size(Prcp,2)
            A(7*size(Prcp,1)+j,Jcp==j)=Rsc;
            Bs(7*size(Prcp,1)+j)=-Rsp(:,j);
        end
        %         [~,indipendent_lineA]=licols(full(A)');
        %         A=sparse(A(indipendent_lineA,:));
        %         Bs=sparse(Bs(indipendent_lineA,:));
        A=sparse(A);
        Bs=sparse(Bs);
        x0=Prcp(:);
        x0=sparse(x0);
        Lambda=2*((A*A')\(Bs-A*x0));
        x=x0+A'*Lambda/2;
        Prcp=reshape(x,size(Prcp,1),size(Prcp,2));
        A=zeros(7*size(Prpc,1)+size(Prpc,2),length(Prcp(:)));
        Bs=zeros(7*size(Prpc,1)+size(Prpc,2),1);
        %        A=zeros(6*size(Prpc,1),length(Prcp(:)));
        %         Bs=zeros(6*size(Prpc,1),1);
        for i=1:size(Prpc,1)
            A((i-1)*7+(1:7),Ipc==i)=[As2;B2(:,3*(i-1)+(1:3))';Dc'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            Bs((i-1)*7+(1:7))=[b;Dp(i,:)'];
        end
        for j=1:size(Prpc,2)
            A(7*size(Prpc,1)+j,Jpc==j)=Rsp;
            Bs(7*size(Prpc,1)+j)=-Rsc(:,j);
        end
        %        [~,indipendent_lineA]=licols(full(A)');
        %         A=sparse(A(indipendent_lineA,:));
        %         Bs=sparse(Bs(indipendent_lineA,:));
        A=sparse(A);
        Bs=sparse(Bs);
        x0=Prpc(:);
        x0=sparse(x0);
        Lambda=2*((A*A')\(Bs-A*x0));
        x=x0+A'*Lambda/2;
        Prpc=reshape(x,size(Prpc,1),size(Prpc,2));
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==5
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        %analytic solution constant stress
        for nc=1:size(interface_c_coord,1)
            Dc(3*(nc-1)+(1:3),1:3)=repmat(interface_c_coord(nc,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1]+10;
        end
        for np=1:size(interface_p_coord,1)
            Dp(3*(np-1)+(1:3),1:3)=repmat(interface_p_coord(np,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1]+10;
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        Rsc=As2*Mc; Rsp=As*Mp;
        [Icp,Jcp]=find(ones(size(Prcp)));
        [Ipc,Jpc]=find(ones(size(Prpc)));
        A=zeros(6*size(Prcp,1)+3*size(Prcp,2),length(Prcp(:)));
        Bs=zeros(6*size(Prcp,1)+3*size(Prcp,2),1);
        %         A=zeros(6*size(Prcp,1),length(Prcp(:)));
        %         Bs=zeros(6*size(Prcp,1),1);
        for i=1:size(Prcp,1)
            A((i-1)*6+(1:6),Icp==i)=[As;B(:,3*(i-1)+(1:3))'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            Bs((i-1)*6+(1:6))=[b];
        end
        for j=1:size(Prcp,2)
            A(6*size(Prcp,1)+(j-1)*3+(1:3),Jcp==j)=Rsc;
            Bs(6*size(Prcp,1)+(j-1)*3+(1:3))=Rsp(:,j);
        end
        %         [~,indipendent_lineA]=licols(full(A)');
        %         A=sparse(A(indipendent_lineA,:));
        %         Bs=sparse(Bs(indipendent_lineA,:));
        A=sparse(A);
        Bs=sparse(Bs);
        x0=Prcp(:);
        x0=sparse(x0);
        Lambda=2*((A*A')\(Bs-A*x0));
        x=x0+A'*Lambda/2;
        Prcp=reshape(x,size(Prcp,1),size(Prcp,2));
        %        A=zeros(9*size(Prpc,1)+3*size(Prpc,2),length(Prcp(:)));
        %         Bs=zeros(9*size(Prpc,1)+3*size(Prpc,2),1);
        A=zeros(6*size(Prpc,1)+3*size(Prpc,2),length(Prcp(:)));
        Bs=zeros(6*size(Prpc,1)+3*size(Prpc,2),1);
        for i=1:size(Prpc,1)
            A((i-1)*6+(1:6),Ipc==i)=[As2;B2(:,3*(i-1)+(1:3))'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            Bs((i-1)*6+(1:6))=[b];
        end
        for j=1:size(Prpc,2)
            A(6*size(Prpc,1)+(j-1)*3+(1:3),Jpc==j)=Rsp;
            Bs(6*size(Prpc,1)+(j-1)*3+(1:3))=Rsc(:,j);
        end
        %        [~,indipendent_lineA]=licols(full(A)');
        %         A=sparse(A(indipendent_lineA,:));
        %         Bs=sparse(Bs(indipendent_lineA,:));
        A=sparse(A);
        Bs=sparse(Bs);
        x0=Prpc(:);
        x0=sparse(x0);
        Lambda=2*((A*A')\(Bs-A*x0));
        x=x0+A'*Lambda/2;
        Prpc=reshape(x,size(Prpc,1),size(Prpc,2));
        Q_cp=Prpc.'; Q_pc=Prcp.';
    end
    if interp_residual==1
        Q_pc=Mp*(Prcp)'/Mc;Q_cp=Mc*(Prpc)'/Mp;
    end
elseif Intergrid_switch==19
    %Mortar corrigé
    %     [GP_coord,El,detJ]=Gauss_point(FACE_int_p,Plate_coord);
    [csi_eta]=identification_Q4function(Cube_coord,FACE_int_c,interface_p_coord);
    [csi_eta2]=identification_Q4function(Plate_coord,FACE_int_p,interface_c_coord);
    [Mortar_GPpc,mortar_detjpc]=Geometric_intersection_quadratic(interface_c_coord,ELintic,interface_p_coord,ELintip,csi_eta,csi_eta2);
    [Mortar_GPcp,mortar_detjcp]=Geometric_intersection_quadratic(interface_p_coord,ELintip,interface_c_coord,ELintic,csi_eta2,csi_eta);
    %     [Mc]=intergrid_mass_matrix_intersection(Mortar_GPcp,mortar_detjcp,FACE_int_c,FACE_int_p,Cube_DOFs_Map);
    [Mc]=intergrid_mass_matrix_intersection(Mortar_GPpc,mortar_detjpc,FACE_int_c,FACE_int_p,Cube_DOFs_Map);
    
    Mc=Mc(int_dofs_c(:,5),int_dofs_c(:,5));
    %     [Mp]=intergrid_mass_matrix_intersection(Mortar_GPpc,mortar_detjpc,FACE_int_p,FACE_int_c,Plate_DOFs_Map);
    [Mp]=intergrid_mass_matrix_intersection(Mortar_GPcp,mortar_detjcp,FACE_int_p,FACE_int_c,Plate_DOFs_Map);
    
    Mp=Mp(int_dofs_p(:,5),int_dofs_p(:,5));
    [Mpc]=intergrid_mass_matrix_special(Mortar_GPcp,Mortar_GPpc,mortar_detjcp,FACE_int_p,FACE_int_c,Cube_DOFs_Map,Plate_DOFs_Map);
    [Mcp]=intergrid_mass_matrix_special(Mortar_GPpc,Mortar_GPcp,mortar_detjpc,FACE_int_c,FACE_int_p,Plate_DOFs_Map,Cube_DOFs_Map);
    Mcp=Mcp(int_dofs_c(:,5),int_dofs_p(:,5));
    Mpc=Mpc(int_dofs_p(:,5),int_dofs_c(:,5));
    Prpc=Mp\Mpc;Prcp=Mc\Mcp;
    Prpc=diag(sum(Mp,2))*Prpc/(diag(sum(Mc,2))); Prcp=diag(sum(Mc,2))*Prcp/(diag(sum(Mp,2)));
    Q_cp=Prpc.'; Q_pc=Prcp.';
    if correzione_momento==1
        %tensore dei bracci.
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        for i=1:size(Prcp,1)
            A=repmat((Prcp(i,:)~=0)',1,3);
            C=(null(full(B(:,3*(i-1)+(1:3)).*A)'))';
            Prcp(i,:)=(Prcp(i,:)*C.')*C;
            Prcp(i,:)=Prcp(i,:)/(sum(Prcp(i,:))*(sum(Prcp(i,:))~=0)+1*(sum(Prcp(i,:))==0));
        end
        for i=1:size(Prpc,1)
            A=repmat((Prpc(i,:)~=0)',1,3);
            C=(null(full(B2(:,3*(i-1)+(1:3)).*A)'))';
            Prpc(i,:)=(Prpc(i,:)*C.')*C;
            Prpc(i,:)=Prpc(i,:)/(sum(Prpc(i,:))*(sum(Prpc(i,:))~=0)+1*(sum(Prpc(i,:))==0));
        end
        %         C=(null((B)'))';
        %         C2=(null((B2)'))';
        %         Prcp=(Prcp*C.')*C;Prpc=(Prpc*C2.')*C2;
        %         Prcp=diag(sum(Prcp,2))\Prcp;
        %         Prpc=diag(sum(Prpc,2))\Prpc;
        
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==2
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        for i=1:size(Prcp,1)
            A=[As;B(:,3*(i-1)+(1:3))'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            x0=Prcp(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prcp(i,:)=x';
        end
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        for i=1:size(Prpc,1)
            A=[As2;B2(:,3*(i-1)+(1:3))'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            x0=Prpc(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prpc(i,:)=x';
        end
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==3
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        %analytic solution constant stress
        for nc=1:size(interface_c_coord,1)
            Dc(nc,1:3)=repmat(interface_c_coord(nc,:),3,1).*[-.29 -.29 1];
        end
        for np=1:size(interface_p_coord,1)
            Dp(np,1:3)=repmat(interface_p_coord(np,:),3,1).*[-.29 -.29 1];
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        for i=1:size(Prcp,1)
            A=[As;B(:,3*(i-1)+(1:3))';Dp'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            b=[b;Dc(i,:)'];
            x0=Prcp(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prcp(i,:)=x';
        end
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        for i=1:size(Prpc,1)
            A=[As2;B2(:,3*(i-1)+(1:3))';Dc'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            b=[b;Dp(i,:)'];
            x0=Prpc(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prpc(i,:)=x';
        end
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==4
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        %analytic solution constant stress
        for nc=1:size(interface_c_coord,1)
            Dc(3*(nc-1)+(1:3),1:3)=repmat(interface_c_coord(nc,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1];
        end
        for np=1:size(interface_p_coord,1)
            Dp(3*(np-1)+(1:3),1:3)=repmat(interface_p_coord(np,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1];
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        Rsc=As2*Mc; Rsp=As*Mp;
        [Icp,Jcp]=find(ones(size(Prcp)));
        [Ipc,Jpc]=find(ones(size(Prpc)));
        A=zeros(9*size(Prcp,1)+3*size(Prcp,2),length(Prcp(:)));
        Bs=zeros(9*size(Prcp,1)+3*size(Prcp,2),1);
        for i=1:size(Prcp,1)
            A((i-1)*9+(1:9),Icp==i)=[As;B(:,3*(i-1)+(1:3))';Dp'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            Bs((i-1)*9+(1:9))=[b;Dc(i,:)'];
        end
        for j=1:size(Prcp,2)
            A(9*size(Prcp,1)+(j-1)*3+(1:3),Jcp==j)=Rsc;
            Bs(9*size(Prcp,1)+(j-1)*3+(1:3))=Rsp(:,j);
        end
        [~,indipendent_lineA]=licols(full(A)');
        A=sparse(A(indipendent_lineA,:));
        Bs=sparse(Bs(indipendent_lineA,:));
        x0=Prcp(:);
        x0=sparse(x0);
        Lambda=2*((A*A')\(Bs-A*x0));
        x=x0+A'*Lambda/2;
        Prcp=reshape(x,size(Prcp,1),size(Prcp,2));
        A=zeros(9*size(Prcp,1)+3*size(Prcp,2),length(Prcp(:)));
        Bs=zeros(9*size(Prcp,1)+3*size(Prcp,2),1);
        for i=1:size(Prpc,1)
            A((i-1)*9+(1:9),Icp==i)=[As2;B2(:,3*(i-1)+(1:3))';Dc'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            Bs((i-1)*9+(1:9))=[b;Dp(i,:)'];
        end
        for j=1:size(Prpc,2)
            A(9*size(Prpc,1)+(j-1)*3+(1:3),Jpc==j)=Rsp;
            Bs(9*size(Prpc,1)+(j-1)*3+(1:3))=Rsc(:,j);
        end
        [~,indipendent_lineA]=licols(full(A)');
        A=sparse(A(indipendent_lineA,:));
        Bs=sparse(Bs(indipendent_lineA,:));
        x0=Prpc(:);
        x0=sparse(x0);
        Lambda=2*((A*A')\(Bs-A*x0));
        x=x0+A'*Lambda/2;
        Prpc=reshape(x,size(Prpc,1),size(Prpc,2));
        Q_cp=Prpc.'; Q_pc=Prcp.';
        
    end
end

% sorting of DOFs for intrerfacing
[cube_only_dof]=setdiff(Cube_DOFs_Map(:,4),int_dofs_c(:,5),'stable');
[plate_only_dof]=setdiff(Plate_DOFs_Map(:,4),int_dofs_p(:,5),'stable');

% [Ic,Jc,k_ijc,Isc,Jsc,DBsc]=Stiffness_assembly(Cube_coord,Ec,CUBE_ZONE_ELEMENT_DOFS);
%
% [Ip,Jp,k_ijp,Isp,Jsp,DBsp]=Stiffness_assembly(Plate_coord,Ep,PLATE_ZONE_ELEMENT_DOFS);

% BC
hold on
clamped_dofs=[3*clamped_nodes-2,3*clamped_nodes-1,3*clamped_nodes];
% clamped_dofs=[3*clamped_nodes];
% clamped_dofs=[clamped_dofs;3*(total_clamped)-2;3*(total_clamped)-1;3*(yzclamped)-1];
plot3(Plate_coord(clamped_nodes,1),Plate_coord(clamped_nodes,2),Plate_coord(clamped_nodes,3),'k^','MarkerFaceColor','k')
% clamped_dofs=[clamped_dofs;3*clamped_nodes(1)-2;3*clamped_nodes(1)-1;3*clamped_nodes(5)-1];
% clamped_dofs=identification_vector_plate(clamped_dofs(:));
activedofs=setdiff((1:(3*length(Xp))'),clamped_dofs);

%find the condensated matrix of the plate
% Kp=sparse(Ip,Jp,k_ijp,3*length(Xp),3*length(Xp)); Kp=Kp+triu(Kp.',1);
% int_dof_plate=int_dofs_p(:,5);
observation_dof_plate=setdiff((1:(3*length(Xp))'),int_dof_plate);
[~,clamped_dofs_id]=intersect(observation_dof_plate,clamped_dofs(:));
Kpoo=Kp(observation_dof_plate,observation_dof_plate);
observation_dof_plate_not_clamp=setdiff(observation_dof_plate,observation_dof_plate(clamped_dofs_id));
Kpoo(clamped_dofs_id,:)=[];
Kpoo(:,clamped_dofs_id)=[];
Kpoi=Kp(observation_dof_plate,int_dof_plate);
Kpoi(clamped_dofs_id,:)=[];
Kpio=Kpoi.';
Kpii=Kp(int_dof_plate,int_dof_plate);
if primal_schwartz==1
    RMp=-Kpoo\Kpoi;
    Kpii_tilder=Kpii+Kpio*RMp;
end


%find the condensated matrix of the cube
% Kc=sparse(Ic,Jc,k_ijc,3*length(Xc),3*length(Xc)); Kc=Kc+triu(Kc.',1);
% int_dof_cube=int_dofs_c(:,5);
observation_dof_cube=setdiff((1:(3*length(Xc))'),int_dof_cube);
Kcoo=Kc(observation_dof_cube,observation_dof_cube);
Kcoi=Kc(observation_dof_cube,int_dof_cube);
Kcio=Kcoi.';
Kcii=Kc(int_dof_cube,int_dof_cube);
if primal_schwartz==1
    RMc=-Kcoo\Kcoi;
    Kcii_tilder=Kcii+Kcio*RMc;
end
%Nodal or constant pressure application

if pressure_switch==1
    
    [ID4]=ismember(FACES_c(:,[2:5]),excited_nodes);
    IDF=ID4(:,1)&ID4(:,2)&ID4(:,3)&ID4(:,4);
    FACE_pressure=FACES_c(IDF,:);
    load_direction=Load_dir;
    [fp]=constant_unit_pressure_load(FACE_pressure,Cube_DOFs_Map,Cube_coord,load_direction);
    % Resultant = 121N on z
    constant_pressure=121/4;
    Fc=constant_pressure*fp;
else
    %excitation on the cube
    excited_dofs=[3*excited_nodes];
    excited_dofs=excited_dofs(:);
    Fc=sparse(excited_dofs,ones(size(excited_dofs)),ones(size(excited_dofs)),3*length(Xc),1);
end
%load partition
Fco=Fc(observation_dof_cube);
Fci=Fc(int_dof_cube);
if primal_schwartz==1
    %load reduction
    Fci_tilder=Fci+RMc.'*Fco;
    u0c=Kcoo\Fco;
end

For_rep=full(reshape(Fc,3,[]))';
hold on
quiver3(Xc(excited_nodes),Yc(excited_nodes),Zc(excited_nodes),For_rep(excited_nodes,1),For_rep(excited_nodes,2),For_rep(excited_nodes,3),'LineWidth',1)

% Displacement in both cases

if (Intergrid_switch~=12&&Intergrid_switch~=13)&&(Intergrid_switch~=14)&&Intergrid_switch~=15&&(Intergrid_switch~=18)
    if cube_master==1
        if primal_schwartz==1
            %when the cube is master:
            Kci_tot=Kcii_tilder+Q_cp*Kpii_tilder*Prpc;
            Fci_tot=Fci_tilder;
            uci=Kci_tot\Fci_tot;
            upi=Prpc*uci;
            uco=u0c+RMc*uci;
            upo=RMp*upi;
        end
        K_tot=[Kcoo,Kcoi,zeros(size(Kcoi,1),size(Kpoo,2))
            Kcio Kcii+Q_cp*Kpii*Prpc Q_cp*Kpio
            zeros(size(Kcoi,1),size(Kpoo,2))' Kpoi*Prpc Kpoo];
        F=[Fco;Fci;zeros(size(Kpoo,2),1)];
        U=K_tot\F;
        uco=U(1:length(Fco));
        uci=U(length(Fco)+(1:length(Fci)));
        upi=Prpc*uci;
        upo=U(length(Fco)+length(Fci)+(1:size(Kpoo,2)));
    else
        if primal_schwartz==1
            %when the plate is master:
            Kpi_tot=Kpii_tilder+Q_pc*Kcii_tilder*Prcp;
            Fpi_tot=Q_pc*Fci_tilder;
            upi=Kpi_tot\Fpi_tot;
            uci=Prcp*upi;
            uco=u0c+RMc*uci;
            upo=RMp*upi;
        end
        K_tot=[Kcoo,Kcoi*Prcp,zeros(size(Kcoi,1),size(Kpoo,2))
            Q_pc*Kcio Kpii+Q_pc*Kcii*Prcp Kpio
            zeros(size(Kcoi,1),size(Kpoo,2))' Kpoi Kpoo];
        F=[Fco;Q_pc*Fci;zeros(size(Kpoo,2),1)];
        U=K_tot\F;
        uco=U(1:length(Fco));
        upi=U(length(Fco)+(1:size(Kpii,1)));
        uci=Prcp*upi;
        upo=U(length(Fco)+size(Kpii,1)+(1:size(Kpoo,2)));
    end
    if nomaster==1
        K=[Kcii_tilder+2*Kcii_tilder*eye(size(Kcii_tilder))+2*Prpc'*Kpii_tilder*Prpc,-2*(Kcii_tilder*Prcp+Prpc'*Kpii_tilder);-2*(Kpii_tilder*Prpc+Prcp'*Kcii_tilder),Kpii_tilder+2*Kpii_tilder*eye(size(Kpii_tilder))+2*Prcp'*Kcii_tilder*Prcp];
        F=[Fci_tilder;zeros(size(Kpii_tilder,1),1)];
        U=K\F;
        uci=U(1:size(Kcii_tilder,1));
        upi=U((size(Kcii_tilder,1)+1):end);
        uco=u0c+RMc*uci;
        upo=RMp*upi;
    end
elseif Intergrid_switch==15
    %normal
    Kcpii=Kcp(:,int_dof_plate);
    Kcpio=Kcp(:,observation_dof_plate);
    Kcpio(:,clamped_dofs_id)=[];
    Kpcii=Kpc(:,int_dof_cube);
    Kpcio=Kpc(:,observation_dof_cube);
    Kcicii=Kcic(:,int_dof_cube);
    Kcicio=Kcic(:,observation_dof_cube);
    Kpipii=Kpip(:,int_dof_plate);
    Kpipio=Kpip(:,observation_dof_plate);
    %tangent
    Ktcpii=Ktcp(:,int_dof_plate);
    Ktcpio=Ktcp(:,observation_dof_plate);
    Ktcpio(:,clamped_dofs_id)=[];
    Ktpcii=Ktpc(:,int_dof_cube);
    Ktpcio=Ktpc(:,observation_dof_cube);
    Ktcicii=Ktcic(:,int_dof_cube);
    Ktcicio=Ktcic(:,observation_dof_cube);
    Ktpipii=Ktpip(:,int_dof_plate);
    Ktpipio=Ktpip(:,observation_dof_plate);
    %bitangent
    Kbcpii=Kbcp(:,int_dof_plate);
    Kbcpio=Kbcp(:,observation_dof_plate);
    Kbcpio(:,clamped_dofs_id)=[];
    Kbpcii=Kbpc(:,int_dof_cube);
    Kbpcio=Kbpc(:,observation_dof_cube);
    Kbcicii=Kbcic(:,int_dof_cube);
    Kbcicio=Kbcic(:,observation_dof_cube);
    Kbpipii=Kbpip(:,int_dof_plate);
    Kbpipio=Kbpip(:,observation_dof_plate);
    %Stiffness assembly
    K_tot=[[Kcoo,Kcoi;Kcio,Kcii],[Kcicio';Kcicii'],zeros(size([Kcoo,Kcoi;Kcio,Kcii],1),size([Kpii Kpio;Kpoi Kpoo],2)),;
        [Kcicio,Kcicii],zeros(size(Kcicio,1),size(Kcicio,1)),[Kcpii,Kcpio];
        zeros(size([Kcoo,Kcoi;Kcio,Kcii],1),size([Kpii Kpio;Kpoi Kpoo],2))',[Kcpii,Kcpio]',[Kpii Kpio;Kpoi Kpoo]];
    j1=size(Kcoo,1)+size(Kcio,1)+size(Kcio,1);
    %
    K_tot=[K_tot(1:j1,1:j1),[zeros(size(Kcoo,1),size(Mcp,2));Mcp;zeros(size(Kcii,1),size(Mcp,2))],K_tot(1:j1,(j1+1):end)
        [[zeros(size(Kcoo,1),size(Mcp,2));Mcp;zeros(size(Kcii,1),size(Mcp,2))]' ,zeros(size(Mcp,2)),[-Mp,zeros(size(Mcp,2),size(Kpoo,1))]]
        K_tot((j1+1):end,1:j1),[-Mp,zeros(size(Mcp,2),size(Kpoo,1))]',K_tot((j1+1):end,(j1+1):end)];
    %
    j2=size(Kcoo,1)+size(Kcio,1);
    K_tot=[[K_tot(1:j2,1:j2),[Ktcicio';Ktcicii'],K_tot(1:j2,(j2+1):end)];
        [[Ktcicio,Ktcicii],zeros(size(Ktcicio,1),size(K_tot,2)+size([Ktcicio';Ktcicii'],2)-size([Ktcicio,Ktcicii],2)-size([Ktcpii,Ktcpio],2)),[Ktcpii,Ktcpio]];
        [K_tot((j2+1):end,1:j2),[zeros(size(K_tot((j2+1):end,1:j2),1)-size([Ktcpii,Ktcpio]',1),size([Ktcpii,Ktcpio]',2));[Ktcpii,Ktcpio]'],K_tot((j2+1):end,(j2+1):end)]];
    %
    K_tot=[[K_tot(1:j2,1:j2),[Kbcicio';Kbcicii'],K_tot(1:j2,(j2+1):end)];
        [[Kbcicio,Kbcicii],zeros(size(Kbcicio,1),size(K_tot,2)+size([Kbcicio';Kbcicii'],2)-size([Kbcicio,Kbcicii],2)-size([Kbcpii,Kbcpio],2)),[Kbcpii,Kbcpio]];
        [K_tot((j2+1):end,1:j2),[zeros(size(K_tot((j2+1):end,1:j2),1)-size([Kbcpii,Kbcpio]',1),size([Kbcpii,Kbcpio]',2));[Kbcpii,Kbcpio]'],K_tot((j2+1):end,(j2+1):end)]];
    %
    F=[Fco;Fci;zeros(size(Kpoo,2)+size(Kpii,2)+size(Kcicio,1),1)];
    F=[F(1:j1);
        zeros(size(Mcp,2),1);
        F((j1+1):end)];
    F=[F(1:j2);
        zeros(size(Ktcicio,1),1);
        F((j2+1):end)];
    F=[F(1:j2);
        zeros(size(Kbcicio,1),1);
        F((j2+1):end)];
    U=K_tot\F;
    uco=U(1:length(Fco));
    uci=U(length(Fco)+(1:size(Kcii,1)));
    lambdab=U(length(Fco)+size(Kcii,1)+(1:size(Kcii,1)));
    lambdat=U(length(Fco)+2*size(Kcii,1)+(1:size(Kcii,1)));
    lambdan=U(length(Fco)+3*size(Kcii,1)+(1:size(Kcii,1)));
    psip=U(length(Fco)+4*size(Kcii,1)+(1:size(Kpii,2)));
    upi=U(length(Fco)+4*size(Kcii,1)+size(Kpii,2)+(1:size(Kpii,2)));
    upo=U(length(Fco)+4*size(Kcii,1)+2*size(Kpii,2)+(1:size(Kpoo,2)));
elseif Intergrid_switch==14
    Kcpii=Kcp(:,int_dof_plate);
    Kcpio=Kcp(:,observation_dof_plate);
    Kcpio(:,clamped_dofs_id)=[];
    Kpcii=Kpc(:,int_dof_cube);
    Kpcio=Kpc(:,observation_dof_cube);
    Kcicii=Kcic(:,int_dof_cube);
    Kcicio=Kcic(:,observation_dof_cube);
    Kpipii=Kpip(:,int_dof_plate);
    Kpipio=Kpip(:,observation_dof_plate);
    K_tot=[[Kcoo,Kcoi;Kcio,Kcii],[Kcicio';Kcicii'],zeros(size([Kcoo,Kcoi;Kcio,Kcii],1),size([Kpii Kpio;Kpoi Kpoo],2)),;
        [Kcicio,Kcicii],zeros(size(Kcicio,1),size(Kcicio,1)),[Kcpii,Kcpio];
        zeros(size([Kcoo,Kcoi;Kcio,Kcii],1),size([Kpii Kpio;Kpoi Kpoo],2))',[Kcpii,Kcpio]',[Kpii Kpio;Kpoi Kpoo]];
    j1=size(Kcoo,1)+size(Kcio,1)+size(Kcio,1);
    K_tot=[K_tot(1:j1,1:j1),[zeros(size(Kcoo,1),size(Mcp,2));Mcp;zeros(size(Kcii,1),size(Mcp,2))],K_tot(1:j1,(j1+1):end)
        [[zeros(size(Kcoo,1),size(Mcp,2));Mcp;zeros(size(Kcii,1),size(Mcp,2))]' ,zeros(size(Mcp,2)),[-Mp,zeros(size(Mcp,2),size(Kpoo,1))]]
        K_tot((j1+1):end,1:j1),[-Mp,zeros(size(Mcp,2),size(Kpoo,1))]',K_tot((j1+1):end,(j1+1):end)];
    F=[Fco;Fci;zeros(size(Kpoo,2)+size(Kpii,2)+size(Kcicio,1),1)];
    F=[F(1:j1);
        zeros(size(Mcp,2),1);
        F((j1+1):end)];
    U=K_tot\F;
    uco=U(1:length(Fco));
    uci=U(length(Fco)+(1:size(Kcii,1)));
    lambdac=U(length(Fco)+size(Kcii,1)+(1:size(Kcii,1)));
    psip=U(length(Fco)+2*size(Kcii,1)+(1:size(Kpii,2)));
    upi=U(length(Fco)+2*size(Kcii,1)+size(Kpii,2)+(1:size(Kpii,2)));
    upo=U(length(Fco)+2*size(Kcii,1)+2*size(Kpii,2)+(1:size(Kpoo,2)));
elseif Intergrid_switch==13
    %     [Vc,Dc]=eig(full(Kcii));
    %     Dc=diag(Dc);
    %     [~,indc]=sort(Dc);
    %     Vc=Vc(:,indc);
    %     Defc=Vc(:,7:end);
    %     [Vp,Dp]=eig(full(Kpii));
    %     Defp=Vp(:,7:end);
    Kcpii=Kcp(:,int_dof_plate);
    Kcpio=Kcp(:,observation_dof_plate);
    Kcpio(:,clamped_dofs_id)=[];
    Kpcii=Kpc(:,int_dof_cube);
    Kpcio=Kpc(:,observation_dof_cube);
    Kcicii=Kcic(:,int_dof_cube);
    Kcicio=Kcic(:,observation_dof_cube);
    Kpipii=Kpip(:,int_dof_plate);
    Kpipio=Kpip(:,observation_dof_plate);
    K_tot=[[Kcoo,Kcoi;Kcio,Kcii],-[Kcicio';Kcicii'],zeros(size([Kcoo,Kcoi;Kcio,Kcii],1),size([Kpii Kpio;Kpoi Kpoo],2));
        -[Kcicio,Kcicii],zeros(size(Kcicio,1),size(Kcicio,1)),[Kcpii,Kcpio];
        zeros(size([Kcoo,Kcoi;Kcio,Kcii],1),size([Kpii Kpio;Kpoi Kpoo],2))',[Kcpii,Kcpio]',[Kpii Kpio;Kpoi Kpoo]];
    F=[Fco;Fci;zeros(size(Kpoo,2)+size(Kpii,2)+size(Kcicio,1),1)];
    Ico=1:size(Fco,1);
    Ici=size(Fco,1)+(1:size(Fci,1));
    Icd=size(Fco,1)+(1:(size(Fci,1)-6));
    Icl1=size(Fco,1)+size(Fci,1)+(1:size(Fci,1));
    Ipi=Icl1(end)+(1:size(Kpii,1));
    Ipo1=Ipi(end)+(1:(size(Kpoo,2)));
    Icl2=Icl1-6;
    Ial=Icl2(end)+(1:6);
    Ipd=Ial(end)+(1:(size(Kpii,2)-6));
    Ipo2=Ipd(end)+(1:(size(Kpoo,2)));
    Tr=zeros(size(F,1),size(F,1)-6);
    Tr(Ico,Ico)=speye(size(Fco,1));
    Tr(Ici,[Icd,Ial])=[Defc,Rigid_c];
    Tr(Icl1,Icl2)=speye(length(Icl1));
    Tr(Ipi,[Ipd,Ial])=[Defp,Rigid_p];
    Tr(Ipo1,Ipo2)=speye(length(Ipo1));
    %     figure
    %     spy(Tr)
    K_tilder_tot=Tr'*K_tot*Tr;
    F_tilder_tot=Tr'*F;
    Utilder=K_tilder_tot\F_tilder_tot;
    U=Tr*Utilder;
    %     U=K_tot\F;
    uco=U(1:length(Fco));
    uci=U(length(Fco)+(1:size(Kcii,1)));
    upi=U(length(Fco)+size(Kcii,1)+size(Kcicio,1)+(1:size(Kpii,2)));
    upo=U(length(Fco)+size(Kcii,1)+size(Kpii,2)+size(Kcicio,1)+(1:size(Kpoo,2)));
elseif Intergrid_switch==18
    if correzione_momento==1
        %tensore dei bracci.
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        for i=1:size(Prcp,1)
            A=repmat((Prcp(i,:)~=0)',1,3);
            C=(null((B(:,3*(i-1)+(1:3)).*A)'))';
            Prcp(i,:)=(Prcp(i,:)*C.')*C;
            Prcp(i,:)=Prcp(i,:)/(sum(Prcp(i,:))*(sum(Prcp(i,:))~=0)+1*(sum(Prcp(i,:))==0));
        end
        for i=1:size(Prpc,1)
            A=repmat((Prpc(i,:)~=0)',1,3);
            C=(null((B2(:,3*(i-1)+(1:3)).*A)'))';
            Prpc(i,:)=(Prpc(i,:)*C.')*C;
            Prpc(i,:)=Prpc(i,:)/(sum(Prpc(i,:))*(sum(Prpc(i,:))~=0)+1*(sum(Prpc(i,:))==0));
        end
        %         C=(null((B)'))';
        %         C2=(null((B2)'))';
        %         Prcp=(Prcp*C.')*C;Prpc=(Prpc*C2.')*C2;
        %         Prcp=diag(sum(Prcp,2))\Prcp;
        %         Prpc=diag(sum(Prpc,2))\Prpc;
        
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==2
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        for i=1:size(Prcp,1)
            A=[As;B(:,3*(i-1)+(1:3))'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            x0=Prcp(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prcp(i,:)=x';
        end
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        for i=1:size(Prpc,1)
            A=[As2;B2(:,3*(i-1)+(1:3))'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            x0=Prpc(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prpc(i,:)=x';
        end
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==3
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        %analytic solution constant stress
        for nc=1:size(interface_c_coord,1)
            Dc(3*(nc-1)+(1:3),1:3)=repmat(interface_c_coord(nc,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1];
        end
        for np=1:size(interface_p_coord,1)
            Dp(3*(np-1)+(1:3),1:3)=repmat(interface_p_coord(np,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1];
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        for i=1:size(Prcp,1)
            A=[As;B(:,3*(i-1)+(1:3))';Dp'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            b=[b;Dc(i,:)'];
            x0=Prcp(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prcp(i,:)=x';
        end
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        for i=1:size(Prpc,1)
            A=[As2;B2(:,3*(i-1)+(1:3))';Dc'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            b=[b;Dp(i,:)'];
            x0=Prpc(i,:)';
            Lambda=2*((A*A')\(b-A*x0));
            x=x0+A'*Lambda/2;
            Prpc(i,:)=x';
        end
        Q_cp=Prpc.'; Q_pc=Prcp.';
    elseif correzione_momento==4
        for nc=1:size(interface_c_coord,1)
            for np=1:size(interface_p_coord,1)
                P2P1=interface_c_coord(nc,:)-interface_p_coord(np,:);
                P1P2=-P2P1;
                for j=1:3
                    for i=1:3
                        di=zeros(1,3);
                        dj=zeros(1,3);
                        di(i)=1;
                        dj(j)=1;
                        B(3*(np-1)+j,3*(3*(nc-1)+i-1)+(1:3))=(cross(P1P2,dj));
                        B2(3*(nc-1)+i,3*(3*(np-1)+j-1)+(1:3))=(cross(P2P1,di));
                    end
                end
            end
        end
        %analytic solution constant stress
        for nc=1:size(interface_c_coord,1)
            Dc(3*(nc-1)+(1:3),1:3)=repmat(interface_c_coord(nc,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1];
        end
        for np=1:size(interface_p_coord,1)
            Dp(3*(np-1)+(1:3),1:3)=repmat(interface_p_coord(np,:),3,1).*[1 -0.29 -0.29;-0.29 1 -.29;-.29 -.29 1];
        end
        As=repmat(eye(3),1,size(Prcp,2)/3);
        As2=repmat(eye(3),1,size(Prpc,2)/3);
        Rsc=As2*Mc; Rsp=As*Mp;
        [Icp,Jcp]=find(ones(size(Prcp)));
        [Ipc,Jpc]=find(ones(size(Prpc)));
        A=zeros(9*size(Prcp,1)+3*size(Prcp,2),length(Prcp(:)));
        Bs=zeros(9*size(Prcp,1)+3*size(Prcp,2),1);
        for i=1:size(Prcp,1)
            A((i-1)*9+(1:9),Icp==i)=[As;B(:,3*(i-1)+(1:3))';Dp'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            Bs((i-1)*9+(1:9))=[b;Dc(i,:)'];
        end
        for j=1:size(Prcp,2)
            A(9*size(Prcp,1)+(j-1)*3+(1:3),Jcp==j)=Rsc;
            Bs(9*size(Prcp,1)+(j-1)*3+(1:3))=Rsp(:,j);
        end
        [~,indipendent_lineA]=licols(full(A)');
        A=sparse(A(indipendent_lineA,:));
        Bs=sparse(Bs(indipendent_lineA,:));
        x0=Prcp(:);
        x0=sparse(x0);
        Lambda=2*((A*A')\(Bs-A*x0));
        x=x0+A'*Lambda/2;
        Prcp=reshape(x,size(Prcp,1),size(Prcp,2));
        A=zeros(9*size(Prcp,1)+3*size(Prcp,2),length(Prcp(:)));
        Bs=zeros(9*size(Prcp,1)+3*size(Prcp,2),1);
        for i=1:size(Prpc,1)
            A((i-1)*9+(1:9),Icp==i)=[As2;B2(:,3*(i-1)+(1:3))';Dc'];
            b=zeros(6,1);
            b(rem(i-1,3)+1)=1;
            Bs((i-1)*9+(1:9))=[b;Dp(i,:)'];
        end
        for j=1:size(Prpc,2)
            A(9*size(Prpc,1)+(j-1)*3+(1:3),Jpc==j)=Rsp;
            Bs(9*size(Prpc,1)+(j-1)*3+(1:3))=Rsc(:,j);
        end
        [~,indipendent_lineA]=licols(full(A)');
        A=sparse(A(indipendent_lineA,:));
        Bs=sparse(Bs(indipendent_lineA,:));
        x0=Prpc(:);
        x0=sparse(x0);
        Lambda=2*((A*A')\(Bs-A*x0));
        x=x0+A'*Lambda/2;
        Prpc=reshape(x,size(Prpc,1),size(Prpc,2));
        Q_cp=Prpc.'; Q_pc=Prcp.';
    end
    Total_DOFs=length(observation_dof_plate_not_clamp)+length(observation_dof_cube)+length(int_dof_plate)+length(int_dof_cube);
    Total_multiplicators=length(int_dof_plate)+length(int_dof_cube);
    Active_dof_cube=size([Kcoo,Kcoi;Kcio,Kcii],1);
    Active_dof_plate=size([Kpoo,Kpoi;Kpio,Kpii],1);
    G1=[-eye(size(Prpc,2)),Prcp;Prpc,-eye(size(Prpc,1))];
    [~,indipendent_line]=licols(full(G1)');
    dipendent_line=setdiff(1:size(G1,1),indipendent_line);
    K_tot=zeros(Total_DOFs+Total_multiplicators);
    K_tot(Active_dof_cube+(1:length(int_dof_cube)),length(observation_dof_cube)+(1:length(int_dof_cube)))=-eye(length(int_dof_cube));
    K_tot(Active_dof_cube+(1:length(int_dof_cube)),Active_dof_cube+Total_multiplicators+(1:length(int_dof_plate)))=Prcp;
    K_tot(Active_dof_cube+Total_multiplicators+(1:length(int_dof_plate)),Active_dof_cube+(1:length(int_dof_cube)))=Prcp';
    K_tot(length(observation_dof_cube)+(1:length(int_dof_cube)),Active_dof_cube+(1:length(int_dof_cube)))=-eye(length(int_dof_cube));
    K_tot(Active_dof_cube+length(int_dof_cube)+(1:length(int_dof_plate)),Active_dof_cube+Total_multiplicators+(1:length(int_dof_plate)))=-eye(length(int_dof_plate));
    K_tot(Active_dof_cube+Total_multiplicators+(1:length(int_dof_plate)),Active_dof_cube+length(int_dof_cube)+(1:length(int_dof_plate)))=-eye(length(int_dof_plate));
    K_tot(Active_dof_cube+length(int_dof_cube)+(1:length(int_dof_plate)),length(observation_dof_cube)+(1:length(int_dof_cube)))=Prpc;
    K_tot(length(observation_dof_cube)+(1:length(int_dof_cube)),Active_dof_cube+length(int_dof_cube)+(1:length(int_dof_plate)))=Prpc';
    K_tot(1:Active_dof_cube,1:Active_dof_cube)=[Kcoo,Kcoi;Kcio,Kcii];
    K_tot((end-Active_dof_plate+1):end,(end-Active_dof_plate+1):end)=[Kpii Kpio;Kpoi Kpoo];
    
    F=[Fco;Fci;zeros(Active_dof_plate+Total_multiplicators,1)];
    K_tot(Active_dof_cube+dipendent_line,:)=[];
    K_tot(:,Active_dof_cube+dipendent_line)=[];
    F(Active_dof_cube+dipendent_line)=[];
    % eliminate innecessary multiplicators and equations
    U=K_tot\F;
    uco=U(1:length(Fco));
    uci=U(length(Fco)+(1:size(Kcii,1)));
    upi=U(Active_dof_cube+Total_multiplicators-length(dipendent_line)+(1:length(int_dof_plate)));
    upo=U(Active_dof_cube+Total_multiplicators+length(int_dof_plate)-length(dipendent_line)+(1:length(observation_dof_plate_not_clamp)));
elseif Intergrid_switch==12
    Kcpii=Kcp(:,int_dof_plate);
    Kcpio=Kcp(:,observation_dof_plate);
    Kcpio(:,clamped_dofs_id)=[];
    Kpcii=Kpc(:,int_dof_cube);
    Kpcio=Kpc(:,observation_dof_cube);
    K_tot=[[Kcoo,Kcoi,
        Kcio,Kcii] ,[zeros(size(Kcoo,1),size(Kcpii,2)+size(Kcpio,2))
        Kcpii Kcpio];
        [Kpcio Kpcii
        zeros(size(Kpoi,1),size(Kpcio,2)+size(Kpcii,2))],[Kpii Kpio;Kpoi Kpoo]];
    F=[Fco;Fci;zeros(size(Kpoo,2)+size(Kpii,2),1)];
    U=K_tot\F;
    %     %filter the rigid body motion
    %     [Vc,Dc]=eigs([Kcoo,Kcoi;Kcio,Kcii],6,'SM');
    % %     [Vp,Dp]=eigs(zeros(size([Kpii Kpio;Kpoi Kpoo])),6,'SM');
    %     Vp=zeros(size([Kpii Kpio;Kpoi Kpoo],1),6);
    %     V=[Vc;Vp];
    %     U=U - sum((ones(size(U))*(U'*V)).*V,2);
    uco=U(1:length(Fco));
    uci=U(length(Fco)+(1:size(Kcii,1)));
    upi=U(length(Fco)+size(Kcii,1)+(1:size(Kpii,2)));
    upo=U(length(Fco)+size(Kcii,1)+size(Kpii,2)+(1:size(Kpoo,2)));
end
Ucube_c([observation_dof_cube(:);int_dof_cube],1)=[uco;uci];

Uplate_c([observation_dof_plate_not_clamp(:);clamped_dofs(:);int_dof_plate],1)=[upo;zeros(size(clamped_dofs(:)));upi];

Grid_motion_cube_c=zeros(length(Xc),3);
Grid_motion_plate_c=zeros(length(Xp),3);
% Average dimension
for nn=1:size(Cube_DOFs_Map,1)
    Grid_motion_cube_c(Cube_DOFs_Map(nn,1),Cube_DOFs_Map(nn,3))=Ucube_c(Cube_DOFs_Map(nn,4));
end
for nn=1:size(Plate_DOFs_Map,1)
    Grid_motion_plate_c(Plate_DOFs_Map(nn,1),Plate_DOFs_Map(nn,3))=Uplate_c(Plate_DOFs_Map(nn,4));
end
Grid_motion=[Grid_motion_cube_c;Grid_motion_plate_c];
% Plot
COORD=[Cube_coord;Plate_coord];
MeanSize = mean(max(COORD)-min(COORD));
% Maximum motion
[MaxGridMotion,I] = max(abs(real(Grid_motion(:))));
% Color displacement
ColorOdsR = sqrt(sum(real(Grid_motion)'.^2))';
% ColorOdsR = Grid_motion(:,3);


Final_coord=zeros(size(COORD));
for k=1:size(COORD,1)
    Final_coord(k,:)=COORD(k,:)+0.25*MeanSize*Grid_motion(k,:)/MaxGridMotion;
end
FACES=[FACES_c;FACES_p+length(Xc)];
figure
hold on; h25_patch = patch('Vertices',Final_coord(1:length(Xc),:),'Faces',FACES_c(:,2:end),'CData',ColorOdsR(1:length(Xc)),'FaceColor','interp','SpecularColorReflectance',0.3);
hold on; h25_patch = patch('Vertices',Final_coord(length(Xc)+(1:length(Xp)),:),'Faces',FACES_p(:,2:end),'CData',ColorOdsR(length(Xc)+(1:length(Xp))),'FaceColor','interp','SpecularColorReflectance',0.3);
colormap jet;
quiver3(0,0,0,4,0,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,4,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,0,4,'LineWidth',3,'Color','k')
text(4,0,0,'x')
text(0,4,0,'y')
text(0,0,4,'z')
view([135 35]);
axis equal
axis off
c=colorbar;
c.Limits=[min(sqrt(sum(real(Grid_motion)'.^2))') max(sqrt(sum(real(Grid_motion)'.^2))')];
caxis([min(sqrt(sum(real(Grid_motion)'.^2))') max(sqrt(sum(real(Grid_motion)'.^2))')])
c.Label.String = 'displacement [mm] ';
c.FontWeight='bold';
%plot of the interface displacement
figure
subplot(2,1,1)
hold on; h25_patch = patch('Vertices',Cube_coord,'Faces',FACE_int_c(:,2:end),'CData',ColorOdsR(1:length(Xc)),'FaceColor','interp','SpecularColorReflectance',0.3);
c=colorbar;
colormap jet;
c.Limits=[min(sqrt(sum(real(Grid_motion)'.^2))') max(sqrt(sum(real(Grid_motion)'.^2))')];
caxis([min(sqrt(sum(real(Grid_motion)'.^2))') max(sqrt(sum(real(Grid_motion)'.^2))')])
c.Label.String = ' displacement cube interface [mm]  ';
c.FontWeight='bold';
axis equal
axis off
%
subplot(2,1,2)
hold on; h25_patch = patch('Vertices',Plate_coord,'Faces',FACE_int_p(:,2:end),'CData',ColorOdsR(length(Xc)+(1:length(Xp))),'FaceColor','interp','SpecularColorReflectance',0.3);
c=colorbar;
colormap jet;
c.Limits=[min(sqrt(sum(real(Grid_motion)'.^2))') max(sqrt(sum(real(Grid_motion)'.^2))')];
caxis([min(sqrt(sum(real(Grid_motion)'.^2))') max(sqrt(sum(real(Grid_motion)'.^2))')])
c.Label.String = ' displacement plate interface [mm] ';
c.FontWeight='bold';
axis equal
axis off
% constraint [MPa]


% use the element extrapolation and simple smoothening
if average_stress==1&&global_least_sq==0
    Ips=zeros(6*size(Ep,1)*size(Ep,2),1);
    Jps=Ips;
    pps=ones(size(Ips));
    for k=1:size(Ep,1)
        for j=1:size(Ep,2)
            Ips((k-1)*size(Ep,2)*6+(j-1)*6+(1:6))=6*(Ep(k,j)-1)+(1:6);
            Jps((k-1)*size(Ep,2)*6+(j-1)*6+(1:6))=(k-1)*size(Ep,2)*6+(j-1)*6+(1:6);
        end
    end
    Ps=sparse(Ips,Jps,pps,6*max(Ep(:)),48*size(Ep,1));
    Ls=diag(sum(Ps,2));
    Ssup=Ls\Ps*sparse(Isp,Jsp,DBsp,48*size(Ep,1),3*max(Ep(:)));
    Snp=Ssup*Uplate_c;
    Snp=reshape(Snp,6,[])';
    s1=Snp(:,1);
    s2=Snp(:,2);
    s3=Snp(:,3);
    s12=Snp(:,4);
    s23=Snp(:,5);
    s31=Snp(:,6);
    Snp=sqrt(0.5*((s1-s2).^2+(s2-s3).^2+(s3-s1).^2+6*(s12.^2+s23.^2+s31.^2)));
    Ips=zeros(6*size(Ec,1)*size(Ec,2),1);
    Jps=Ips;
    pps=ones(size(Ips));
    for k=1:size(Ec,1)
        for j=1:size(Ec,2)
            Ips((k-1)*size(Ec,2)*6+(j-1)*6+(1:6))=6*(Ec(k,j)-1)+(1:6);
            Jps((k-1)*size(Ec,2)*6+(j-1)*6+(1:6))=(k-1)*size(Ec,2)*6+(j-1)*6+(1:6);
        end
    end
    Ps=sparse(Ips,Jps,pps,6*max(Ec(:)),48*size(Ec,1));
    Ls=diag(sum(Ps,2));
    Ssuc=Ls\Ps*sparse(Isc,Jsc,DBsc,48*size(Ec,1),3*max(Ec(:)));
    Snc=Ssuc*Ucube_c;
    Snc=reshape(Snc,6,[])';
    s1=Snc(:,1);
    s2=Snc(:,2);
    s3=Snc(:,3);
    s12=Snc(:,4);
    s23=Snc(:,5);
    s31=Snc(:,6);
    Snc=sqrt(0.5*((s1-s2).^2+(s2-s3).^2+(s3-s1).^2+6*(s12.^2+s23.^2+s31.^2)));
elseif average_stress==1&&global_least_sq==1
    Sgpp=sparse(Isp,Jsp,DBsp,48*size(Ep,1),3*max(Ep(:)))*Uplate_c;
    Nnp=sparse(Inp,Jnp,Nnnp,size(Sgpp,1),6*max(Ep(:)));
    stress_not_int_p=[6*(not_int_p(:)-1)+1,6*(not_int_p(:)-1)+2,6*(not_int_p(:)-1)+3,6*(not_int_p(:)-1)+4,6*(not_int_p(:)-1)+5,6*(not_int_p(:)-1)+6];
    stress_not_int_p=reshape(stress_not_int_p',[],1);
    stress_int_p=[6*(interface_p(:)-1)+1,6*(interface_p(:)-1)+2,6*(interface_p(:)-1)+3,6*(interface_p(:)-1)+4,6*(interface_p(:)-1)+5,6*(interface_p(:)-1)+6];
    stress_int_p=reshape(stress_int_p',[],1);
    stress_not_int_c=[6*(not_int_c(:)-1)+1,6*(not_int_c(:)-1)+2,6*(not_int_c(:)-1)+3,6*(not_int_c(:)-1)+4,6*(not_int_c(:)-1)+5,6*(not_int_c(:)-1)+6];
    stress_not_int_c=reshape(stress_not_int_c',[],1);
    stress_int_c=[6*(interface_c(:)-1)+1,6*(interface_c(:)-1)+2,6*(interface_c(:)-1)+3,6*(interface_c(:)-1)+4,6*(interface_c(:)-1)+5,6*(interface_c(:)-1)+6];
    stress_int_c=reshape(stress_int_c',[],1);
%     Snp=(Nnp'*Nnp)\(Nnp'*Sgpp);
   
    Sgpc=sparse(Isc,Jsc,DBsc,48*size(Ec,1),3*max(Ec(:)))*Ucube_c;
    Nnc=sparse(Inc,Jnc,Nnnc,size(Sgpc,1),6*max(Ec(:)));
    Nn1=Nnc(:,stress_not_int_c);
    Nng1=Nnc(:,stress_int_c);
    Nn2=Nnp(:,stress_not_int_p);
    Nng2=Nnp(:,stress_int_p);
    Nn_11g22g=[Nn1,Nng1,zeros(size(Nng1,1),size(Nnp,2));zeros(size(Nng2,1),size(Nnc,2)),Nng2,Nn2];
    I1=speye(6*length(not_int_c));
    I2=speye(6*length(not_int_p));
    Ig1=speye(6*length(interface_c));
    Ig2=speye(6*length(interface_p));
    [Ipicp,Jpicp,Prcp_noden]=find(Pr_node_cp);
    Istresscp=[6*(Ipicp-1)+1,6*(Ipicp-1)+2,6*(Ipicp-1)+3,6*(Ipicp-1)+4,6*(Ipicp-1)+5,6*(Ipicp-1)+6];
    Istresscp=reshape(Istresscp',[],1);
    Jstresscp=[6*(Jpicp-1)+1,6*(Jpicp-1)+2,6*(Jpicp-1)+3,6*(Jpicp-1)+4,6*(Jpicp-1)+5,6*(Jpicp-1)+6];
    Jstresscp=reshape(Jstresscp',[],1);
    Prcp_noden=reshape(repmat(Prcp_noden,1,6)',[],1);
    Pr_cp_stress=sparse(Istresscp,Jstresscp,Prcp_noden,6*size(Pr_node_cp,1),6*size(Pr_node_cp,2));
    [Ipipc,Jpipc,Prpc_noden]=find(Pr_node_pc);
    Istresspc=[6*(Ipipc-1)+1,6*(Ipipc-1)+2,6*(Ipipc-1)+3,6*(Ipipc-1)+4,6*(Ipipc-1)+5,6*(Ipipc-1)+6];
    Istresspc=reshape(Istresspc',[],1);
    Jstresspc=[6*(Jpipc-1)+1,6*(Jpipc-1)+2,6*(Jpipc-1)+3,6*(Jpipc-1)+4,6*(Jpipc-1)+5,6*(Jpipc-1)+6];
    Jstresspc=reshape(Jstresspc',[],1);
    Prpc_noden=reshape(repmat(Prpc_noden,1,6)',[],1);
    Pr_pc_stress=sparse(Istresspc,Jstresspc,Prpc_noden,6*size(Pr_node_pc,1),6*size(Pr_node_pc,2));
    Tcm=[I1,zeros(size(I1,1),size(Ig1,2)),zeros(size(I1,1),size(I2,2))
        zeros(size(Ig1,1),size(I1,2)),Ig1,zeros(size(Ig1,1),size(I2,2))
        zeros(size(Ig2,1),size(I1,2)),Pr_pc_stress,zeros(size(Ig2,1),size(I2,2))
        zeros(size(I2,1),size(I1,2)),zeros(size(I2,1),size(Ig1,2)),I2];
    Tpm=[I1,zeros(size(I1,1),size(Ig2,2)),zeros(size(I1,1),size(I2,2))
        zeros(size(Ig1,1),size(I1,2)),Pr_cp_stress,zeros(size(Ig1,1),size(I2,2))
        zeros(size(Ig2,1),size(I1,2)),Ig2,zeros(size(Ig2,1),size(I2,2))
        zeros(size(I2,1),size(I1,2)),zeros(size(I2,1),size(Ig2,2)),I2];
    Nn_11g2=Nn_11g22g*Tcm;Nn_12g2=Nn_11g22g*Tpm;
    if cube_master==1
        Snodal=Tcm*((Nn_11g2'*Nn_11g2)\(Nn_11g2'*[Sgpc;Sgpp]));
    else
        Snodal=Tpm*((Nn_12g2'*Nn_12g2)\(Nn_12g2'*[Sgpc;Sgpp]));
    end
%     Snc=zeros(size(Sgpp));
    Snc([stress_not_int_c;stress_int_c])=Snodal(1:6*length(Zc));
    Snp([stress_int_p;stress_not_int_p])=Snodal(6*length(Zc)+(1:6*length(Zp)));
    %     Snc=(Nnc'*Nnc)\(Nnc'*Sgpc);
    Snp=reshape(Snp,6,[])';
    s1=Snp(:,1);
    s2=Snp(:,2);
    s3=Snp(:,3);
    s12=Snp(:,4);
    s23=Snp(:,5);
    s31=Snp(:,6);
    Snp=sqrt(0.5*((s1-s2).^2+(s2-s3).^2+(s3-s1).^2+6*(s12.^2+s23.^2+s31.^2)));
    Snc=reshape(Snc,6,[])';
    s1=Snc(:,1);
    s2=Snc(:,2);
    s3=Snc(:,3);
    s12=Snc(:,4);
    s23=Snc(:,5);
    s31=Snc(:,6);
    Snc=sqrt(0.5*((s1-s2).^2+(s2-s3).^2+(s3-s1).^2+6*(s12.^2+s23.^2+s31.^2)));
else
    SUp=sparse(Isp,Jsp,DBsp,6*size(Ep,1),3*length(Xp));
    Selp=SUp*Uplate_c;
    con_sel=3;
    In=zeros(8*size(Ep,1),1);
    Jel=In;
    cont=In;
    for k=1:size(Ep,1)
        In(8*(k-1)+(1:8))=Ep(k,:);
        Jel(8*(k-1)+(1:8))=6*(k-1)+con_sel;
        cont(8*(k-1)+(1:8))=1;
    end
    N_ELp=sparse(In,Jel,cont,length(Xp),6*size(Ep,1));
    Snp=(N_ELp*Selp)./(N_ELp*ones(size(Selp)));
    SUc=sparse(Isc,Jsc,DBsc,6*size(Ec,1),3*length(Xc));
    Selc=SUc*Ucube_c;
    In=zeros(8*size(Ec,1),1);
    Jel=In;
    cont=In;
    for k=1:size(Ec,1)
        In(8*(k-1)+(1:8))=Ec(k,:);
        Jel(8*(k-1)+(1:8))=6*(k-1)+con_sel;
        cont(8*(k-1)+(1:8))=1;
    end
    N_ELc=sparse(In,Jel,cont,length(Xc),6*size(Ec,1));
    Snc=N_ELc*Selc./(N_ELc*ones(size(Selc)));
    if Vhon_mises==1
        Snc=[];
        Snp=[];
        for con_sel=1:6
            SUp=sparse(Isp,Jsp,DBsp,6*size(Ep,1),3*length(Xp));
            Selp=SUp*Uplate_c;
            %     con_sel=3;
            In=zeros(8*size(Ep,1),1);
            Jel=In;
            cont=In;
            for k=1:size(Ep,1)
                In(8*(k-1)+(1:8))=Ep(k,:);
                Jel(8*(k-1)+(1:8))=6*(k-1)+con_sel;
                cont(8*(k-1)+(1:8))=1;
            end
            N_ELp=sparse(In,Jel,cont,length(Xp),6*size(Ep,1));
            Snp=[Snp,(N_ELp*Selp)./(N_ELp*ones(size(Selp)))];
            SUc=sparse(Isc,Jsc,DBsc,6*size(Ec,1),3*length(Xc));
            Selc=SUc*Ucube_c;
            In=zeros(8*size(Ec,1),1);
            Jel=In;
            cont=In;
            for k=1:size(Ec,1)
                In(8*(k-1)+(1:8))=Ec(k,:);
                Jel(8*(k-1)+(1:8))=6*(k-1)+con_sel;
                cont(8*(k-1)+(1:8))=1;
            end
            N_ELc=sparse(In,Jel,cont,length(Xc),6*size(Ec,1));
            Snc=[Snc,N_ELc*Selc./(N_ELc*ones(size(Selc)))];
        end
        s1=Snc(:,1);
        s2=Snc(:,2);
        s3=Snc(:,3);
        s12=Snc(:,4);
        s23=Snc(:,5);
        s31=Snc(:,6);
        Snc=sqrt(0.5*((s1-s2).^2+(s2-s3).^2+(s3-s1).^2+6*(s12.^2+s23.^2+s31.^2)));
        s1=Snp(:,1);
        s2=Snp(:,2);
        s3=Snp(:,3);
        s12=Snp(:,4);
        s23=Snp(:,5);
        s31=Snp(:,6);
        Snp=sqrt(0.5*((s1-s2).^2+(s2-s3).^2+(s3-s1).^2+6*(s12.^2+s23.^2+s31.^2)));
    end
end
% Snpn=Snp;Sncn=Snc;
% Snpn(interface_p)=1/2*(Snp(interface_p)+Pr_node_pc*Snc(interface_c));%averaging the interface stress
% Sncn(interface_c)=1/2*(Snc(interface_c)+Pr_node_cp*Snp(interface_p));
% Snp=Snpn;Snc=Sncn;
ColorOdsR = [Snc;Snp];
% ColorOdsR = ColorOdsR - min(ColorOdsR); ColorOdsR = ColorOdsR/max(ColorOdsR);
figure
hold on; h25_patch = patch('Vertices',Final_coord(1:length(Xc),:),'Faces',FACES_c(:,2:end),'CData',ColorOdsR(1:length(Xc)),'FaceColor','interp','SpecularColorReflectance',0.3);
hold on; h25_patch = patch('Vertices',Final_coord(length(Xc)+(1:length(Xp)),:),'Faces',FACES_p(:,2:end),'CData',ColorOdsR(length(Xc)+(1:length(Xp))),'FaceColor','interp','SpecularColorReflectance',0.3);
colormap jet;
quiver3(0,0,0,4,0,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,4,0,'LineWidth',3,'Color','k')
quiver3(0,0,0,0,0,4,'LineWidth',3,'Color','k')
text(4,0,0,'x','FontSize',18,'FontWeight','Bold')
text(0,4,0,'y','FontSize',18,'FontWeight','Bold')
text(0,0,4,'z','FontSize',18,'FontWeight','Bold')
view([135 35]);
c=colorbar;
% caxis([min([Snc;Snp]) max(([Snc;Snp]))])
% c.Limits=[min([Snc;Snp]) max(([Snc;Snp]))];
caxis([0 700])
c.Limits=[0 700];
c.Label.String = '\sigma_{VM} [MPa] ';
c.FontWeight='bold';
axis equal
axis off
if correzione_momento==0
print(['is',num2str(Intergrid_switch)],'-dpng')
elseif correzione_momento==2
    print(['is',num2str(Intergrid_switch),'c'],'-dpng')
end
%plot of the interface Stress
figure
subplot(2,1,1)
hold on; h25_patch = patch('Vertices',Cube_coord,'Faces',FACE_int_c(:,2:end),'CData',ColorOdsR(1:length(Xc)),'FaceColor','interp','SpecularColorReflectance',0.3);
c=colorbar;
colormap jet;
caxis([min([Snc;Snp]) max(([Snc;Snp]))])
c.Limits=[min([Snc;Snp]) max(([Snc;Snp]))];
c.Label.String = '\sigma_{VM} cube interface [MPa]  ';
c.FontWeight='bold';
axis equal
axis off
%
subplot(2,1,2)
hold on; h25_patch = patch('Vertices',Plate_coord,'Faces',FACE_int_p(:,2:end),'CData',ColorOdsR(length(Xc)+(1:length(Xp))),'FaceColor','interp','SpecularColorReflectance',0.3);
c=colorbar;
colormap jet;
caxis([min([Snc;Snp]) max(([Snc;Snp]))])
c.Limits=[min([Snc;Snp]) max(([Snc;Snp]))];
c.Label.String = '\sigma_{VM} plate interface [MPa] ';
c.FontWeight='bold';
axis equal
axis off
%verify the BCs
Rp=Kp*Uplate_c;
Rc=Kc*Ucube_c;
Rplb=Rp(clamped_dofs);
%resultant on downside of the Plate
Force=sum(Rplb)
%check of interface
rc1=Rc(int_dof_cube);
rc2=Rp(int_dof_plate);
RC1=zeros(length(int_dof_cube)/3,3);
for k=1:length(int_dof_cube)/3
    RC1(k,:)=[rc1(3*(k-1)+1),rc1(3*(k-1)+2),rc1(3*(k-1)+3)];
end
RC2=zeros(length(int_dof_plate)/3,3);
for k=1:length(int_dof_plate)/3
    RC2(k,:)=[rc2(3*(k-1)+1),rc2(3*(k-1)+2),rc2(3*(k-1)+3)];
end
% Residual on both interface sides [Rx,Ry,RZ]
%cube side
Fc1=sum(RC1)
%plate side
Fc2=sum(RC2)
%momentum cube side
Mc1=sum(cross(interface_c_coord,RC1,2))
%momentum plate side
Mc2=sum(cross(interface_p_coord,RC2,2))
% compliance in 2 ways:
% total compliance (1)+(2)
comp1=Rp.'*Uplate_c+Rc.'*Ucube_c
% interface 1 compliance
cint1=rc1.'*Ucube_c(int_dof_cube)
% interface 2 compliance
cint2=rc2.'*Uplate_c(int_dof_plate)
% interface 1 kinetic specific
ksint1=Ucube_c(int_dof_cube).'*Mc*Ucube_c(int_dof_cube)
% interface 2 kinetic specific
ksint2=Uplate_c(int_dof_plate).'*Mp*Uplate_c(int_dof_plate)
% max displacement
umax=max(sqrt(sum(real(Grid_motion)'.^2))')
% max S33
S33max= max(([Snc;Snp]))
% min S33
S33min= min(([Snc;Snp]))
% Displacement Error estimation
Derr=(rms((Uplate_c(int_dof_plate)-Prpc*Ucube_c(int_dof_cube))))/rms(Ucube_c(int_dof_cube))*100+(rms((Ucube_c(int_dof_cube)-Prcp*Uplate_c(int_dof_plate))))/rms(Uplate_c(int_dof_plate))*100
% % Stress Error estimation
Serr=rms(Snp(interface_p)-Pr_node_pc*Snc(interface_c))/rms(Snp(interface_p))*100+rms(Snc(interface_c)-Pr_node_cp*Snp(interface_p))/rms(Snc(interface_c))*100
%plot of surfacique force
p1=Mc\rc1;
p2=Mp\rc2;
Perr=rms((p2+Prpc*p1))/rms(p2)*100+rms((p1+Prcp*p2))/rms(p1)*100
p1=reshape(p1,3,[])';
p2=reshape(p2,3,[])';
p1mag=sqrt(sum(real(p1.^2)'))';
p2mag=sqrt(sum(real(p2.^2)'))';
ColorOdsR = zeros((length(Xc)+length(Xp)),1);
ColorOdsR([interface_c;length(Xc)+interface_p])=[p1mag;p2mag];
figure
subplot(2,1,1)
hold on; h25_patch = patch('Vertices',Cube_coord,'Faces',FACE_int_c(:,2:end),'CData',ColorOdsR(1:length(Xc)),'FaceColor','interp','SpecularColorReflectance',0.3);
c=colorbar;
colormap jet;
caxis([min([p1mag;p2mag]) max(([p1mag;p2mag]))])
c.Limits=[min([p1mag;p2mag]) max(([p1mag;p2mag]))];
c.Label.String = 'surface tension magnitude [MPa]  ';
c.FontWeight='bold';
axis equal
axis off
%
subplot(2,1,2)
hold on; h25_patch = patch('Vertices',Plate_coord,'Faces',FACE_int_p(:,2:end),'CData',ColorOdsR(length(Xc)+(1:length(Xp))),'FaceColor','interp','SpecularColorReflectance',0.3);
c=colorbar;
colormap jet;
caxis([min([p1mag;p2mag]) max(([p1mag;p2mag]))])
c.Limits=[min([p1mag;p2mag]) max(([p1mag;p2mag]))];
c.Label.String = 'surface tension magnitude [MPa] ';
c.FontWeight='bold';
axis equal
axis off
%plot of absolute residual
R1=reshape(rc1,3,[])';
R2=reshape(rc2,3,[])';
R1mag=sqrt(sum(real(R1.^2)'))';
R2mag=sqrt(sum(real(R2.^2)'))';
ColorOdsR = zeros((length(Xc)+length(Xp)),1);
ColorOdsR([interface_c;length(Xc)+interface_p])=[R1mag;R2mag];
figure
subplot(2,1,1)
hold on; h25_patch = patch('Vertices',Cube_coord,'Faces',FACE_int_c(:,2:end),'CData',ColorOdsR(1:length(Xc)),'FaceColor','interp','SpecularColorReflectance',0.3);
c=colorbar;
colormap jet;
caxis([min([R1mag;R2mag]) max(([R1mag;R2mag]))])
c.Limits=[min([R1mag;R2mag]) max(([R1mag;R2mag]))];
c.Label.String = 'Residual [N]  ';
c.FontWeight='bold';
axis equal
axis off
%
subplot(2,1,2)
hold on; h25_patch = patch('Vertices',Plate_coord,'Faces',FACE_int_p(:,2:end),'CData',ColorOdsR(length(Xc)+(1:length(Xp))),'FaceColor','interp','SpecularColorReflectance',0.3);
c=colorbar;
colormap jet;
caxis([min([R1mag;R2mag]) max(([R1mag;R2mag]))])
c.Limits=[min([R1mag;R2mag]) max(([R1mag;R2mag]))];
c.Label.String = 'residual [N] ';
c.FontWeight='bold';
axis equal
axis off
