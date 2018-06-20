clear
%open the file
ABAQUS_DAT='test.dat';
ABAQUS_DAT_FilePath='/Users/qiuqiuli/Desktop/Abaqus student project';
disp(['> Reading reference ABAQUS [' ABAQUS_DAT '] job in ' ABAQUS_DAT])
fid = fopen([ABAQUS_DAT_FilePath '/' ABAQUS_DAT], 'rt');
% Reads the lines of the file, on after the other
test1=' ss';
NUM_Force_combitation=0;
% store the max stress and their position information
global MAXI_STRESS_COORD ;
MAXI_STRESS_COORD=zeros(36,4);
global COORD_STRESS_100_0
global MAXI_STRESS_POSITION_LIST
MAXI_STRESS_POSITION_LIST=zeros(6000,2);
COORD_STRESS_0_N1000=zeros(10000,3);
COORD_STRESS_100_0=zeros(10000,3);
COORD_STRESS_200_0=zeros(10000,3);
COORD_STRESS_350_0=zeros(10000,3);
COORD_STRESS_575_0=zeros(10000,3);
COORD_STRESS_912_0=zeros(10000,3);
COORD_STRESS_1000_0=zeros(10000,3);
COORD_STRESS_1000_100=zeros(10000,3);
COORD_STRESS_1000_200=zeros(10000,3);
COORD_STRESS_1000_350=zeros(10000,3);
COORD_STRESS_1000_575=zeros(10000,3);
COORD_STRESS_1000_912=zeros(10000,3);
COORD_STRESS_1000_1000=zeros(10000,3);
COORD_STRESS_900_1000=zeros(10000,3);
COORD_STRESS_800_1000=zeros(10000,3);
COORD_STRESS_650_1000=zeros(10000,3);
COORD_STRESS_425_1000=zeros(10000,3);
COORD_STRESS_87_1000=zeros(10000,3);
COORD_STRESS_0_1000=zeros(10000,3);
COORD_STRESS_N100_900=zeros(10000,3);
COORD_STRESS_N200_800=zeros(10000,3);
COORD_STRESS_N350_650=zeros(10000,3);
COORD_STRESS_N575_425=zeros(10000,3);
COORD_STRESS_N912_87=zeros(10000,3);
COORD_STRESS_N1000_0=zeros(10000,3);
COORD_STRESS_N1000_N100=zeros(10000,3);
COORD_STRESS_N1000_N200=zeros(10000,3);
COORD_STRESS_N1000_N350=zeros(10000,3);
COORD_STRESS_N1000_N575=zeros(10000,3);
COORD_STRESS_N1000_N912=zeros(10000,3);
COORD_STRESS_N1000_N1000=zeros(10000,3);
COORD_STRESS_N900_N1000=zeros(10000,3);
COORD_STRESS_N800_N1000=zeros(10000,3);
COORD_STRESS_N650_N1000=zeros(10000,3);
COORD_STRESS_N425_N1000=zeros(10000,3);
COORD_STRESS_N87_N1000=zeros(10000,3);
% save them all in a 3D matrix
global STRESS_COORD;
STRESS_COORD=zeros(10000,3,36) ;

while ~feof(fid) 
    % Reads the lines of the file, on after the other
    Line = [fgetl(fid) blanks(100)];
    % TF1=100, TF2=0
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.100    ,  TOTAL TIME COMPLETED        0.100');
       total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_100_0(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[100,0, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    % TF1=200, TF2=0
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.200    ,  TOTAL TIME COMPLETED        0.200');
       total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_200_0(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[200,0, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
     % TF1=350, TF2=0
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.350    ,  TOTAL TIME COMPLETED        0.350');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
            
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_350_0(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[350,0, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=575, TF2=0
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.575    ,  TOTAL TIME COMPLETED        0.575');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
            
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_575_0(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
       % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[575,0, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=912.5, TF2=0
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.913    ,  TOTAL TIME COMPLETED        0.913');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_912_0(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[912.5,0, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
    % TF1=1000, TF2=0
    if strcmp(Line(1:72),' STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         1.00');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)]; 
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_1000_0(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[1000,0, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=1000, TF2=100
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.100    ,  TOTAL TIME COMPLETED         1.10');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_1000_100(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[1000,100, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=1000, TF2=200
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.200    ,  TOTAL TIME COMPLETED         1.20');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_1000_200(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[1000,200, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=1000, TF2=350
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.350    ,  TOTAL TIME COMPLETED         1.35');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_1000_350(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[1000,350, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=1000, TF2=575
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.575    ,  TOTAL TIME COMPLETED         1.58');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_1000_575(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[1000,575, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
    % TF1=1000, TF2=912.5
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.913    ,  TOTAL TIME COMPLETED         1.91');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_1000_912(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[1000,912.5, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=1000, TF2=1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         2.00');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_1000_1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[1000,1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
        
     % TF1=900, TF2=1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.100    ,  TOTAL TIME COMPLETED         2.10');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_900_1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
        NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[900,1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
    % TF1=800, TF2=1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.200    ,  TOTAL TIME COMPLETED         2.20');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_800_1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[800,1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=650, TF2=1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.350    ,  TOTAL TIME COMPLETED         2.35');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_650_1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[650,1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
        
     % TF1=425, TF2=1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.575    ,  TOTAL TIME COMPLETED         2.58');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_425_1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[425,1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=87.5, TF2=1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.913    ,  TOTAL TIME COMPLETED         2.91');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_87_1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[87.5,1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=0, TF2=1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         3.00');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_0_1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[0,1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-100, TF2=900
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.100    ,  TOTAL TIME COMPLETED         3.10');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N100_900(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-100,900, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-200, TF2=800
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.200    ,  TOTAL TIME COMPLETED         3.20');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N200_800(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-200,800, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-350, TF2=650
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.350    ,  TOTAL TIME COMPLETED         3.35');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N350_650(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-350,650, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
    % TF1=-575, TF2=425
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.575    ,  TOTAL TIME COMPLETED         3.58');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N575_425(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-575,425, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-912.5, TF2=87.5
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.913    ,  TOTAL TIME COMPLETED         3.91');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N912_87(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-912.5,87.5, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
   % TF1=-1000, TF2=0
    if strcmp(Line(1:72),' STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         4.00');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N1000_0(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-1000,0, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-1000, TF2=-100
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.100    ,  TOTAL TIME COMPLETED         4.10');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N1000_N100(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-1000,-100, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-1000, TF2=-200
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.200    ,  TOTAL TIME COMPLETED         4.20');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N1000_N200(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
        % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-1000,-200, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-1000, TF2=-350
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.350    ,  TOTAL TIME COMPLETED         4.35');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N1000_N350(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-1000,-350, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-1000, TF2=-575
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.575    ,  TOTAL TIME COMPLETED         4.57');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N1000_N575(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-1000,-575, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-1000, TF2=-912.5
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.913    ,  TOTAL TIME COMPLETED         4.91');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N1000_N912(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-1000,-912.5, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-1000, TF2=-1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         5.00');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N1000_N1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-1000,-1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-900, TF2=-1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.100    ,  TOTAL TIME COMPLETED         5.10');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N900_N1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-900,-1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-800, TF2=-1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.200    ,  TOTAL TIME COMPLETED         5.20');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N800_N1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-800,-1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
    % TF1=-650, TF2=-1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.350    ,  TOTAL TIME COMPLETED         5.35');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N650_N1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-650,-1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
    % TF1=-425, TF2=-1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.575    ,  TOTAL TIME COMPLETED         5.57');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N425_N1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-425,-1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
     % TF1=-87.5, TF2=-1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED       0.913    ,  TOTAL TIME COMPLETED         5.91');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_N87_N1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[-87.5,-1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
    
    % TF1=0, TF2=-1000
    if strcmp(Line(1:72),' STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         6.00');
        total_nodes=0;
        for i=1:15
            Line = [fgetl(fid) blanks(100)];
        end
        while strcmp(Line(16),'1');
            total_nodes=total_nodes+1;
            COORD_STRESS_0_N1000(total_nodes,:)=[str2double(Line(24:35)),str2double(Line(37:43)),str2double(Line(48:59))];
            Line = [fgetl(fid) blanks(100)];
        end
         % note the maxi stress and position
         NUM_Force_combitation=NUM_Force_combitation+1;
        Line = [fgetl(fid) blanks(100)];
        if strcmp(Line(1:8),' MAXIMUM');
           MAXI_STRESS_value=str2double(Line(23:30));
           Line = [fgetl(fid) blanks(100)];
           MAXI_STRESS_position=str2double(Line(26:30));
           MAXI_STRESS_COORD(NUM_Force_combitation,:)=[0,-1000, MAXI_STRESS_value,MAXI_STRESS_position];
        end
    end
end

% save the node number- coordinate information
MAXI_STRESS_POSITION_LIST(:,1:2)=COORD_STRESS_100_0(1:6000,2:3);

% save every stress distribution
STRESS_COORD(:,:,1)=COORD_STRESS_100_0;
STRESS_COORD(:,:,2)=COORD_STRESS_200_0;
STRESS_COORD(:,:,3)=COORD_STRESS_350_0;
STRESS_COORD(:,:,4)=COORD_STRESS_575_0;
STRESS_COORD(:,:,5)=COORD_STRESS_912_0;
STRESS_COORD(:,:,6)=COORD_STRESS_1000_0;
STRESS_COORD(:,:,7)=COORD_STRESS_1000_100;
STRESS_COORD(:,:,8)=COORD_STRESS_1000_200;
STRESS_COORD(:,:,9)=COORD_STRESS_1000_350;
STRESS_COORD(:,:,10)=COORD_STRESS_1000_575;
STRESS_COORD(:,:,11)=COORD_STRESS_1000_912;
STRESS_COORD(:,:,12)=COORD_STRESS_1000_1000;
STRESS_COORD(:,:,13)=COORD_STRESS_900_1000;
STRESS_COORD(:,:,14)=COORD_STRESS_800_1000;
STRESS_COORD(:,:,15)=COORD_STRESS_650_1000;
STRESS_COORD(:,:,16)=COORD_STRESS_425_1000;
STRESS_COORD(:,:,17)=COORD_STRESS_87_1000;
STRESS_COORD(:,:,18)=COORD_STRESS_0_1000;
STRESS_COORD(:,:,19)=COORD_STRESS_N100_900;
STRESS_COORD(:,:,20)=COORD_STRESS_N200_800;
STRESS_COORD(:,:,21)=COORD_STRESS_N350_650;
STRESS_COORD(:,:,22)=COORD_STRESS_N575_425;
STRESS_COORD(:,:,23)=COORD_STRESS_N912_87;
STRESS_COORD(:,:,24)=COORD_STRESS_N1000_0;
STRESS_COORD(:,:,25)=COORD_STRESS_N1000_N100;
STRESS_COORD(:,:,26)=COORD_STRESS_N1000_N200;
STRESS_COORD(:,:,27)=COORD_STRESS_N1000_N350;
STRESS_COORD(:,:,28)=COORD_STRESS_N1000_N575;
STRESS_COORD(:,:,29)=COORD_STRESS_N1000_N912;
STRESS_COORD(:,:,30)=COORD_STRESS_N1000_N1000;
STRESS_COORD(:,:,31)=COORD_STRESS_N900_N1000;
STRESS_COORD(:,:,32)=COORD_STRESS_N800_N1000;
STRESS_COORD(:,:,33)=COORD_STRESS_N650_N1000;
STRESS_COORD(:,:,34)=COORD_STRESS_N425_N1000;
STRESS_COORD(:,:,35)=COORD_STRESS_N87_N1000;
STRESS_COORD(:,:,36)=COORD_STRESS_0_N1000;