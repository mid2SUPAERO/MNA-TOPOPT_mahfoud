%Gamma matrix
%Substructure containing ns Matrix U->tip clearance stage (s)
%Gamma(s) is a rectangular matrix(N_(s) x length(U))
stage(1).label='FAN';
stage(2).label='IPC1';
stage(3).label='IPC2';
stage(4).label='IPC3';
stage(5).label='HPC1';
stage(6).label='HPC2';
stage(7).label='HPC3';
stage(8).label='HPC4';
stage(9).label='HPC5';
stage(10).label='HPC6';
stage(11).label='HPC7';
stage(12).label='HPC8';
stage(13).label='HPC9';
stage(14).label='HPC10';
stage(15).label='HPT1';
stage(16).label='HPT2';
stage(17).label='IPT1';
stage(18).label='IPT2';
stage(19).label='IPT3';
stage(20).label='IPT4';
load('substructure_results')
load('lin_pert10.dat.mat')
load('Blade_heights')
for s=1:20
    G_numb=size(SUBCASE(1).POINTS,2);
    Id_center=zeros(1,G_numb);
    Id_casing=Id_center;
    for k=1:G_numb
        Id_center(k)=strcmp(SUBCASE(1).POINTS(k).label,[stage(s).label,' CENTER']);
        Id_casing(k)=strcmp(SUBCASE(1).POINTS(k).label,[stage(s).label,' CASING']);
    end
    Index_CENTER=find(Id_center);
    Index_CASING=find(Id_casing);
    Center_coordinates=zeros(3,1);
    for d=1:3
        Center_coordinates(d,1)=SUBCASE(1).POINTS(Index_CENTER).coord(d).Mod;
    end
    circ_coordinates=zeros(3,length(Index_CASING));
    GAMMA=zeros(length(Index_CASING),size(Output_DOFs_Map,1));
    for k=1:length(Index_CASING)
        for d=1:3
            circ_coordinates(d,k)=SUBCASE(1).POINTS(Index_CASING(k)).coord(d).Mod;
        end
        
    end
    radius_vect=circ_coordinates-Center_coordinates*ones(1,length(Index_CASING));
    radius_norm=sqrt((diag(radius_vect'*radius_vect)));
    radius_versor=radius_vect./(ones(3,1)*radius_norm.');
    radius_versor=radius_versor(2:3,:);
    RV(s).mat=radius_versor;
    RV(s).ROTOR_DOF_ID=find((Output_DOFs_Map(:,1)==SUBCASE(1).POINTS(Index_CENTER).num));
    for k=1:length(Index_CASING)
        RV(s).STATOR_DOF_ID(k,:)=find(Output_DOFs_Map(:,1)==SUBCASE(1).POINTS(Index_CASING(k)).num);
        GAMMA(k,RV(s).STATOR_DOF_ID(k,1))=-radius_versor(1,k);
        GAMMA(k,RV(s).STATOR_DOF_ID(k,2))=-radius_versor(2,k);
        GAMMA(k,RV(s).ROTOR_DOF_ID(1))=-radius_versor(1,k);
        GAMMA(k,RV(s).ROTOR_DOF_ID(2))=-radius_versor(2,k);
    end
    Gamma(s).mat=GAMMA;
    Gamma(s).N=length(Index_CASING);
    Gamma(s).lambda=1/Blade_height(s)*100; %performance sensitivity
end
fichier='gamma_structure';
eval([ 'save ' fichier ' Gamma']);
disp([ 'LE FICHIER ' fichier '.mat A ETE CREE']);
clear
