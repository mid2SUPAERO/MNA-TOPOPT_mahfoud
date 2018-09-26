function [new_COORD,new_element,NODE_SET_new,FACES_new,Pu]=mesh_refinement3D(COORD,ELEMENT,NODE_SET,FACES)
% get a finer mesh from the first one by 3D bisection
Last_node_id=size(COORD,1);
new_COORD=[COORD;zeros(19*size(ELEMENT,1),3)];
new_element=zeros(8*size(ELEMENT,1),8);
NODE_SET_new.cardinality= NODE_SET.cardinality;
NODE_SET_new.NID=NODE_SET.NID;
No=size(COORD,1);
FACES_new=zeros(36*size(ELEMENT,1),4);
% Pu=zeros(size(new_COORD,1),size(COORD,1));
% Pu(1:size(COORD,1),1:size(COORD,1))=eye(size(COORD,1));
NE=12*size(ELEMENT,1);
NF=6*size(ELEMENT,1);
Ntt=8*size(ELEMENT,1)+4*NF+2*NE+No; %Number of non zero element in Pu
ipu=zeros(Ntt,1);jpu=ipu;ppu=ipu;
ipu(1:No)=1:No;
jpu(1:No)=1:No;
ppu(1:No)=1;
internal_external=zeros(19,8);
Edges=[1 2;2 3;3,4;4 1;1 5;2 6;3 7;4 8;5 6;6 7;7 8;8 5];
Faces=[1 2 3 4;1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;5 6 7 8];
center=[1 2 3 4 5 6 7 8];
for l=1:size(Edges,1)
    internal_external(l,Edges(l,:))=0.5;
end
for l=1:size(Faces,1)
    internal_external(size(Edges,1)+l,Faces(l,:))=0.25;
end
internal_external(size(Edges,1)+size(Faces,1)+1,center)=0.125;
[Ip,Jp,Pp]=find(internal_external);
Np=length(Ip);
for k=1:size(ELEMENT,1)
    External_vertex=ELEMENT(k,:);
    External_vertex_coord=COORD(External_vertex,:);
    
    
%     internal_coord=[(External_vertex_coord(1,:)+External_vertex_coord(2,:))/2;(External_vertex_coord(2,:)+External_vertex_coord(3,:))/2;(External_vertex_coord(3,:)+External_vertex_coord(4,:))/2;(External_vertex_coord(1,:)+External_vertex_coord(4,:))/2;(External_vertex_coord(1,:)+External_vertex_coord(2,:)+External_vertex_coord(3,:)+External_vertex_coord(4,:))/4];
    internal_coord=internal_external*External_vertex_coord;
    ipu(No+Np*(k-1)+(1:Np))=Last_node_id+Ip;
    jpu(No+Np*(k-1)+(1:Np))=External_vertex(Jp);
    ppu(No+Np*(k-1)+(1:Np))=Pp;
%     Pu(Last_node_id+1,[External_vertex(1),External_vertex(2)])=0.5;
%     Pu(Last_node_id+2,[External_vertex(2),External_vertex(3)])=0.5;
%     Pu(Last_node_id+3,[External_vertex(3),External_vertex(4)])=0.5;
%     Pu(Last_node_id+4,[External_vertex(1),External_vertex(4)])=0.5;
%     Pu(Last_node_id+5,External_vertex([1 2 3 4]))=0.25;
    new_COORD(Last_node_id+(1:size(internal_external,1)),:)=internal_coord;
    N1_27id=[External_vertex,Last_node_id+(1:size(internal_external,1))]';
    im1=ismember(External_vertex(1),NODE_SET.NID);
    im2=ismember(External_vertex(2),NODE_SET.NID);
    im3=ismember(External_vertex(3),NODE_SET.NID);
    im4=ismember(External_vertex(4),NODE_SET.NID);
    im5=ismember(External_vertex(5),NODE_SET.NID);
    im6=ismember(External_vertex(6),NODE_SET.NID);
    im7=ismember(External_vertex(7),NODE_SET.NID);
    im8=ismember(External_vertex(8),NODE_SET.NID);
    if im1&&im2&&im3&&im4
        NODE_SET_new.cardinality= NODE_SET_new.cardinality+5;
        NODE_SET_new.NID=[NODE_SET_new.NID;N1_27id([9 10 11 12 21]')];
    end
    if im1&&im2&&im6&&im5
        NODE_SET_new.cardinality= NODE_SET_new.cardinality+5;
        NODE_SET_new.NID=[NODE_SET_new.NID;N1_27id([9 14 17 13 22]')];
    end
    if im2&&im3&&im7&&im6
        NODE_SET_new.cardinality= NODE_SET_new.cardinality+5;
        NODE_SET_new.NID=[NODE_SET_new.NID;N1_27id([10 15 18 14 23]')];
    end
    if im3&&im4&&im8&&im7
        NODE_SET_new.cardinality= NODE_SET_new.cardinality+5;
        NODE_SET_new.NID=[NODE_SET_new.NID;N1_27id([11 16 19 15 24]')];
    end
    if im4&&im1&&im8&&im5
        NODE_SET_new.cardinality= NODE_SET_new.cardinality+5;
        NODE_SET_new.NID=[NODE_SET_new.NID;N1_27id([12 13 20 16 25]')];
    end
    if im5&&im6&&im7&&im8
        NODE_SET_new.cardinality= NODE_SET_new.cardinality+5;
        NODE_SET_new.NID=[NODE_SET_new.NID;N1_27id([17 18 19 20 26]')];
    end
    
    new_element(8*(k-1)+(1:8),:)=N1_27id([ 1  9 21 12 13 22 27 25;
        9  2 10 21 22 14 23 27;
        10  3 11 21 23 15 24 27;
        12 21 11  4 25 27 24 16;
        13 22 27 25  5 17 26 20;
        22 14 23 27 17  6 18 26;
        23 15 24 27 18  7 19 26;
        25 27 24 16 20 26 19  8]);
    fc=[ 1  9 22 13;
        9  2 14 22;
        22 14  6 17;
        13 22 17  5;
        12 21 27 25;
        21 10 23 27
        25 27 26 20
        27 23 18 26
        4 11 24 16
        11  3 15 24
        16 24 19  8
        24 15  7 19
        4 12 25 16
        12  1 13 25
        16 25 20  8
        25 13  5 20
        11 21 27 24
        21  9 22 27
        27 22 17 26
        24 27 26 19
        3 10 23 15
        10  2 14 23
        15 23 18  7
        23 14  6 18
        1  9 21 12
        9  2 10 21
        21 10  3 11
        12 21 11  4
        13 22 27 25
        22 14 23 27
        27 23 15 24
        25 27 24 16
        5 17 26 20
        17  6 18 26
        26 18  7 19
        20 26 19  8];
    FACES_new(36*(k-1)+(1:36),:)=N1_27id(fc);       
    Last_node_id=size(COORD,1)+size(internal_external,1)*k;
% %     Check the connectivity
%     cor=[-1 -1 -1
%         1 -1 -1
%         1 1 -1
%         -1 1 -1
%         -1 -1 1
%         1 -1 1
%         1 1 1
%         -1 1 1
%         0 -1 -1
%         1 0 -1
%         0 1 -1
%         -1 0 -1
%         -1 -1 0
%         1 -1 0
%         1 1 0
%         -1 1 0
%         0 -1 1
%         1 0 1
%         0 1 1
%         -1 0 1
%         0 0 -1
%         0 -1 0
%         1 0 0
%         0 1 0
%         -1 0 0
%         0 0 1
%         0 0 0];
%     
%     scatter3(cor(:,1),cor(:,2),cor(:,3),'ko','fill')
%     for k=1:27
%         text(cor(k,1),cor(k,2),cor(k,3),[' ',num2str(k)])
%     end
%     patch('Vertices',cor,'Faces',fc,'FaceColor','b','FaceAlpha',0.1)
%     
end
FACES_new=unique(FACES_new,'rows','stable');


Pu=sparse(ipu,jpu,ppu,size(new_COORD,1),size(COORD,1));
Pu=diag((sum(Pu,2)))\Pu;
old_new_COORD=new_COORD;
[new_COORD,~,IC] = unique(old_new_COORD,'rows','stable');
% %oldnew_coord=new_COORD(IC,:);
%new_COORD=oldnew_coord(IA,:);
new_element=IC(new_element);
NODE_SET_new.NID=IC(NODE_SET_new.NID);
NODE_SET_new.NID=unique(NODE_SET_new.NID);
NODE_SET_new.cardinality=length(NODE_SET_new.NID);
FACES_new=IC(FACES_new);
Pr_node_cp=Pu;
% Pu=zeros(2*size(new_COORD,1),2*size(COORD,1));
[I,J,PP]=find(Pr_node_cp);
I=IC(I);
% J=IC(J);
% [I,ij]=unique(IC(I),'stable');
% PP=PP(ij);
% % I=IJ(:,1);
% J=J(ij);
In=zeros(3*length(I),1);
Jn=In;
PPn=In;
In(1:3:end-2)=3*(I-1)+1;
In(2:3:end-1)=3*(I-1)+2;
In(3:3:end)=3*(I-1)+3;
Jn(1:3:end-2)=3*(J-1)+1;
Jn(2:3:end-1)=3*(J-1)+2;
Jn(3:3:end)=3*(J-1)+3;
PPn(1:3:end-2)=PP;
PPn(2:3:end-1)=PP;
PPn(3:3:end)=PP;
Pu=sparse(In,Jn,PPn,3*size(new_COORD,1),3*size(COORD,1));
Pu=diag((sum(Pu,2)))\Pu;
% h = figure(4); set(h,'Color',[1 1 1]);
% hold on; h24_patch = patch('Vertices',new_COORD(:,[1,3]),'Faces',new_element,'FaceVertexCData',[1 1 1],'FaceColor','flat'); axis equal; axis off; drawnow;hold off
