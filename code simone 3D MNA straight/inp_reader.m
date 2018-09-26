clear
ABAQUS_INP='sub110.inp';
ABAQUS_INP_FilePath='C:\Users\s.coniglio\Desktop\Sustructure_topopt -4';
disp(['> Reading reference ABAQUS [' ABAQUS_INP '] job in ' ABAQUS_INP])
fid = fopen([ABAQUS_INP_FilePath '/' ABAQUS_INP], 'rt');
flag = 0;
% Reads the lines of the file, on after the other
Line = [fgetl(fid) blanks(length('*Element, type=S4R'))];
% Spot the end of job
flagBeginBulk = strcmp(Line(1:8),'*Heading');
flagEnddata = strcmp(Line(1:9),'*End Step');
BlockSize = 100;
count = 0; Bound = 0;
GID=zeros(10000000,1);
COORD=zeros(10000000,3);
NumbNSET=0;
total_nodes=0;
total_elements=0;
ELEMENT=[];
while ((feof(fid) == 0)&&(flag == 0))
    % Reads the lines of the file, on after the other
    Line = [fgetl(fid) blanks(length('*Element, type=S4R'))];
    count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
    if (Bound~=LowBound)
        Bound = LowBound;
        display(['       [' num2str(Bound) ' lines]'])
    end
    if strcmp(Line(1:5),'*Node');
        Line = [fgetl(fid) blanks(length('*Element, type=S4R'))];
        count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
        if (Bound~=LowBound)
            Bound = LowBound;
            display(['       [' num2str(Bound) ' lines]'])
        end
        while ~strcmp(Line(1),'*')||strcmp(Line(1:5),'*Node')
            if strcmp(Line(1:5),'*Node');
                
                Line = [fgetl(fid) blanks(length('*Element, type=S4R'))];
                count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
                if (Bound~=LowBound)
                    Bound = LowBound;
                    display(['       [' num2str(Bound) ' lines]'])
                end
            end
            total_nodes=total_nodes+1;
            GID(total_nodes)=str2double(Line(1:8));
            COORD(total_nodes,:)=[str2double(Line(10:22)),str2double(Line(24:36)),str2double(Line(38:50))];
            Line = [fgetl(fid) blanks(length('*Element, type=S4R'))];
            count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
            if (Bound~=LowBound)
                Bound = LowBound;
                display(['       [' num2str(Bound) ' lines]'])
            end
        end
    end
    if strcmp(Line(1:5),'*Nset');
        NumbNSET=NumbNSET+1;
        set_name=deblank(Line(13:end));
        cardinality=0;
        NID=fscanf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d');
        count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
        if (Bound~=LowBound)
            Bound = LowBound;
            display(['       [' num2str(Bound) ' lines]'])
        end
        cardinality=length(NID);
        %         while rem(cardinality,16)==0
        %             fid1=fid;
        %             Line = [fgetl(fid1) blanks(8)];
        %             if strcmp(Line(1),'*')
        %                 if ~strcmp(Line(1:5),'*Nset')
        %                 break
        %                 else
        %                   Line = [fgetl(fid) blanks(8)];
        %                 count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
        %                 if (Bound~=LowBound)
        %                     Bound = LowBound;
        %                     display(['       [' num2str(Bound) ' lines]'])
        %                 end
        %                 end
        %             end
        %             clear fid1
        %             NL=fscanf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d');
        %             count = count + 1; LowBound = floor(count/BlockSize)*BlockSize;
        %             if (Bound~=LowBound)
        %                 Bound = LowBound;
        %                 display(['       [' num2str(Bound) ' lines]'])
        %             end
        %             if length(NL)<1
        %                 break
        %             end
        %             NID=[NID;NL];
        %             cardinality=length(NID);
        %             NumbNSET=NumbNSET+1;
        %         end
        NODE_SET(NumbNSET).name=set_name;
        NODE_SET(NumbNSET).NID=NID;
        NODE_SET(NumbNSET).cardinality=cardinality;
    end
    if strcmp(Line(1:length('*Element, type=S4R')),'*Element, type=S4R')
        new_el=fscanf(fid,'%d,%d,%d,%d,%d');
        new_el=reshape(new_el,5,[]).';
        ELEMENT=[ELEMENT;new_el];
        total_elements=size(ELEMENT,1);
    end
    Line=[Line blanks(length('*Elset, elset=interface_'))];
    if strcmp(Line(1:length('*Elset, elset=interface')),'*Elset, elset=interface')
        if strcmp(Line(1:length('*Elset, elset=interface_')),'*Elset, elset=interface_')
        else
        elID=fscanf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d');
        end
    end
end
[~,~,el_id]=intersect(elID,ELEMENT(:,1));
ELint=ELEMENT(el_id,:);
fclose(fid);
fichier='connectivity_interface';
eval([ 'save ' fichier ' ELint ']);
disp([ 'LE FICHIER ' fichier '.mat A ETE CREE']);
clear
