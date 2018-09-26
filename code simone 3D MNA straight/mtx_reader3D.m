clear
ABAQUS_INP='sub3D.mtx';
ABAQUS_INP_FilePath='C:\Users\s.coniglio\Desktop\Sustructure_topopt -4';
disp(['> Reading reference ABAQUS [' ABAQUS_INP '] job in ' ABAQUS_INP])
fid = fopen([ABAQUS_INP_FilePath '/' ABAQUS_INP], 'rt');
flag = 0;
% Reads the lines of the file, on after the other
Line = [fgetl(fid) blanks(8)];
% Spot the end of job
flagBeginBulk = strcmp(Line(1:8),'*Heading');
flagEnddata = strcmp(Line(1:9),'*End Step');
BlockSize = 10000;
count = 0; Bound = 0;
n_subcase=0;
load('st3D.dat.mat')
while ((feof(fid) == 0)&&(flag == 0))
    % Reads the lines of the file, on after the other
    Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
    count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
    if (Bound~=LowBound)
        Bound = LowBound;
        display(['       [' num2str(Bound) ' lines]'])
    end
    if strcmp(Line(1:13),'*USER ELEMENT');
        Number_of_interface_nodes=str2double(Line(23:32)); %this time there are also lagrangian multipliers of RBE3 that have to be added
        Interface_GID=zeros(Number_of_interface_nodes,1);
        Interface_local_GID=zeros(Number_of_interface_nodes,1);
        Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
        count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
        if (Bound~=LowBound)
            Bound = LowBound;
            display(['       [' num2str(Bound) ' lines]'])
        end
        Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
        count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
        if (Bound~=LowBound)
            Bound = LowBound;
            display(['       [' num2str(Bound) ' lines]'])
        end
        for k=1:10
            Interface_GID(k)=str2double(Line(4+12*(k-1):14+12*(k-1)));
            Interface_local_GID(k)=k;
        end
        cardinality=length(Interface_GID);
        line_number=ceil(Number_of_interface_nodes/10);
        for l=1:(line_number-1)
            Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
            count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
            if (Bound~=LowBound)
                Bound = LowBound;
                display(['       [' num2str(Bound) ' lines]'])
            end
            el=10*(k+10<Number_of_interface_nodes)+rem(Number_of_interface_nodes,10)*(k+10>=Number_of_interface_nodes);
            for  m=1:el
                Interface_GID(m+k)=str2double(Line(4+12*(m-1):14+12*(m-1)));
                Interface_local_GID(m+k)=m+k;
            end
            k=k+el;
        end
        %recognize the ABAQUS INTERNAL DOF (lagrangian multipliers of RBE3)
        Number_of_interface_nodes=sum(Interface_GID>0);
        %DOFs_map 1 2 3
        Retaine_DOFs_Number=3*Number_of_interface_nodes+6*sum(Interface_GID==0);
        DOFs_MAP=zeros(Retaine_DOFs_Number,4);
        for k=1:Number_of_interface_nodes
            DOFs_MAP((k-1)*3+(1:3),:)=[Interface_GID(k)*ones(3,1),Interface_local_GID(k)*ones(3,1),[1;2;3],(k-1)*3+[1;2;3]];
        end
        for k=(Number_of_interface_nodes+1):(Number_of_interface_nodes+sum(Interface_GID==0))
            DOFs_MAP(Number_of_interface_nodes*3+6*(k-Number_of_interface_nodes-1)+(1:6),:)=[Interface_GID(k)*ones(6,1),Interface_local_GID(k)*ones(6,1),[1;2;3;4;5;6],Number_of_interface_nodes*3+(k-Number_of_interface_nodes-1)*6+[1;2;3;4;5;6]];
        end
        while ~strcmp(Line(1:length('*MATRIX,TYPE=STIFFNESS')),'*MATRIX,TYPE=STIFFNESS')
            Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
            count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
            if (Bound~=LowBound)
                Bound = LowBound;
                display(['       [' num2str(Bound) ' lines]'])
            end
        end
        disp('Stiffness Matrix reading')
        K_interface=zeros(Retaine_DOFs_Number);
        Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
        count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
        if (Bound~=LowBound)
            Bound = LowBound;
            display(['       [' num2str(Bound) ' lines]'])
        end
        for l=1:Retaine_DOFs_Number
            line_number=ceil(l/4);
            c=0;
            for s=1:line_number;
                el=4*(c+4<=l)+rem(l,4)*(c+4>l);
                for m=1:el
                    K_interface(l,m+c)=str2double(Line(1+22*(m-1):22+22*(m-1)));
                end
                Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
                count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
                if (Bound~=LowBound)
                    Bound = LowBound;
                    display(['       [' num2str(Bound) ' lines]'])
                end
                c=c+el;
            end
        end
        K_interface=K_interface+triu(K_interface.',1);
    end
    if strcmp(Line(1:length('** SUBSTRUCTURE LOAD')),'** SUBSTRUCTURE LOAD')
        disp('Substructure load vectors reading')
        while strcmp(Line(1:length('** SUBSTRUCTURE LOAD')),'** SUBSTRUCTURE LOAD')
            p=0;
            n_subcase=n_subcase+1;
            LOAD_CASE_VECTOR(n_subcase).title=deblank(Line(47:end));
            for k=1:2
                Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
                count = count + 1; LowBound = floor(count/10000)*10000;
                if (Bound~=LowBound)
                    Bound = LowBound;
                    display(['       [' num2str(Bound) ' lines]'])
                end
            end
            
            while length(deblank(Line))>2&&~strcmp(Line(1:length('** SUBSTRUCTURE LOAD')),'** SUBSTRUCTURE LOAD')
                p=p+1;
                LOAD_CASE_VECTOR(n_subcase).GID(p)=str2double(Line(8:13));
                LOAD_CASE_VECTOR(n_subcase).DOFID(p)=str2double(Line(15:20));
                LOAD_CASE_VECTOR(n_subcase).MOD(p)=str2double(Line(22:34));
                Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
                count = count + 1; LowBound = floor(count/10000)*10000;
                if (Bound~=LowBound)
                    Bound = LowBound;
                    display(['       [' num2str(Bound) ' lines]'])
                end
            end
        end
        LOAD_MATRIX=zeros(Retaine_DOFs_Number,n_subcase);
        for lc=1:n_subcase
            for k=1:length(LOAD_CASE_VECTOR(lc).MOD)
                LOAD_MATRIX(LOAD_CASE_VECTOR(lc).GID(k)==DOFs_MAP(:,2)&LOAD_CASE_VECTOR(lc).DOFID(k)==DOFs_MAP(:,3),lc)=LOAD_CASE_VECTOR(lc).MOD(k);
            end
        end
    end
    
    if strcmp(Line(1:length('** RECOVERY MATRIX')),'** RECOVERY MATRIX')
        Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
        count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
        if (Bound~=LowBound)
            Bound = LowBound;
            display(['       [' num2str(Bound) ' lines]'])
        end
        Number_of_tip_nodes=str2double(deblank(Line(32:end)));
        post_GID=zeros(Number_of_tip_nodes,1);
        post_local_GID=zeros(Number_of_tip_nodes,1);
        
        Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
        count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
        if (Bound~=LowBound)
            Bound = LowBound;
            display(['       [' num2str(Bound) ' lines]'])
        end
        Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
        count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
        if (Bound~=LowBound)
            Bound = LowBound;
            display(['       [' num2str(Bound) ' lines]'])
        end
        for k=1:10
            post_GID(k)=str2double(Line(4+12*(k-1):14+12*(k-1)));
            post_local_GID(k)=k;
        end
        cardinality=length(post_GID);
        line_number=ceil(Number_of_tip_nodes/10);
        for l=1:(line_number-1)
            Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
            count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
            if (Bound~=LowBound)
                Bound = LowBound;
                display(['       [' num2str(Bound) ' lines]'])
            end
            el=10*(k+10<Number_of_tip_nodes)+rem(Number_of_tip_nodes,10)*(k+10>=Number_of_tip_nodes);
            for  m=1:el
                post_GID(m+k)=str2double(Line(4+12*(m-1):14+12*(m-1)));
                post_local_GID(m+k)=m+k;
            end
            k=k+el;
        end
        post_DOFs_MAP=zeros(6*Number_of_tip_nodes,4);
        for k=1:Number_of_tip_nodes
            post_DOFs_MAP((k-1)*6+(1:6),:)=[post_GID(k)*ones(6,1),post_local_GID(k)*ones(6,1),[1;2;3;4;5;6],(k-1)*6+[1;2;3;4;5;6]];
        end
        
        disp('Recovery matrix column selection')
        %list of needed DOFs for tip_clearance:
        %1 for casing nodes
        %2,3 for center nodes
        
        
        G_numb=size(SUBCASE(1).POINTS,2);
        Id_center=zeros(G_numb,1);
        Id_casing=Id_center;
        Id_FAN_center=Id_center;
        for k=1:G_numb
            Id_center(k)=strcmp(SUBCASE(1).POINTS(k).label(end-length(' CENTER')+1:end),' CENTER');
            Id_casing(k)=strcmp(SUBCASE(1).POINTS(k).label(end-length(' CASING')+1:end),' CASING');
            Id_FAN_center(k)=strcmp(SUBCASE(1).POINTS(k).label(1:6),'FAN CE');
        end
        Index_FAN_CENTER=find(Id_FAN_center);
        Index_CENTER=find(Id_center);
        Index_CASING=find(Id_casing);
        CENTER_GID=Index_CENTER;
        CASING_GID=Index_CASING;
        for m=1:length(Index_CENTER)
            CENTER_GID(m)=SUBCASE(1).POINTS(Index_CENTER(m)).num;
        end
        for m=1:length(Index_CASING)
            CASING_GID(m)=SUBCASE(1).POINTS(Index_CASING(m)).num;
        end
        FAN_CENTER_GID=SUBCASE(1).POINTS(Index_FAN_CENTER).num;
        Nuber_of_needed_DOFs=2*length(CASING_GID)+2*length(CENTER_GID);
        Output_DOFs_Map=zeros(Nuber_of_needed_DOFs,4);
        for k=1:size(CASING_GID,1)
            % CT is not possible for condensation to be done in the
            % Retained DOFs regions.
            %
            Output_DOFs_Map(k,:)=post_DOFs_MAP(post_DOFs_MAP(:,1)==CASING_GID(k)&post_DOFs_MAP(:,3)==2,:);
            Output_DOFs_Map(k+size(CASING_GID,1),:)=post_DOFs_MAP(post_DOFs_MAP(:,1)==CASING_GID(k)&post_DOFs_MAP(:,3)==3,:);
        end
        for k=2*size(CASING_GID,1)+1:2*size(CASING_GID,1)+length(CENTER_GID)
            Output_DOFs_Map(k,:)=post_DOFs_MAP(post_DOFs_MAP(:,1)==CENTER_GID(k-2*size(CASING_GID,1))&post_DOFs_MAP(:,3)==2,:);
            Output_DOFs_Map(k+length(CENTER_GID),:)=post_DOFs_MAP(post_DOFs_MAP(:,1)==CENTER_GID(k-2*size(CASING_GID,1))&post_DOFs_MAP(:,3)==3,:);
        end
        [~,sorted_id]=sort(Output_DOFs_Map(:,4));
        Output_DOFs_Map=Output_DOFs_Map(sorted_id,:);
        % add a vector Output this time just for FAN casing Axial
        % displacement
        FAN_center_DOF_MAP=post_DOFs_MAP(post_DOFs_MAP(:,1)==FAN_CENTER_GID&post_DOFs_MAP(:,3)==1,:);
        Recovery_Matrix=zeros(Nuber_of_needed_DOFs,Retaine_DOFs_Number);
        FAN_Center_Recovery_vector=zeros(1,Retaine_DOFs_Number);
    end
    if strcmp(Line(1:length('** SUBSTRUCTURE RECOVERY VECTOR')),'** SUBSTRUCTURE RECOVERY VECTOR')
        %starting a new line of recovery matrix
        %
        % Retained DOF
        Retained_DOF=str2double(deblank(Line(72:end)));
        disp(['Recovery matrix reading retained DOF ',deblank(Line(72:end)),' of ',num2str(Retaine_DOFs_Number)])
       line_number=ceil(6*Number_of_tip_nodes/4);
        c=0;
        for l=1:line_number
            Line = [fgetl(fid) blanks(length('** SUBSTRUCTURE RECOVERY VECTOR'))];
            count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
            if (Bound~=LowBound)
                Bound = LowBound;
                display(['       [' num2str(Bound) ' lines]'])
            end
            el=4*(c+4<=6*Number_of_tip_nodes)+rem(6*Number_of_tip_nodes,4)*(c+4>6*Number_of_tip_nodes);
            for  m=1:el
                if any(c+m==Output_DOFs_Map(:,4))
                    Recovery_Matrix(c+m==Output_DOFs_Map(:,4),Retained_DOF)=str2double(Line(5+22*(m-1):25+22*(m-1)));
                elseif any(c+m==FAN_center_DOF_MAP(:,4))
                    FAN_Center_Recovery_vector(Retained_DOF)=str2double(Line(5+22*(m-1):25+22*(m-1)));
                end
            end
            c=c+el;
        end
    end
end
fichier='substructure_results';
eval([ 'save ' fichier ' Recovery_Matrix Output_DOFs_Map LOAD_MATRIX K_interface DOFs_MAP SUBCASE FAN_Center_Recovery_vector FAN_center_DOF_MAP']);
disp([ 'LE FICHIER ' fichier '.mat A ETE CREE']);
fclose(fid);
clear

