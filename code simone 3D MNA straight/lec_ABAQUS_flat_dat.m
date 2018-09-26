function [SUBCASE] = lec_ABAQUS_flat_dat()
tic
%===========================================================================%
%  LEC_SOL101 Fonction qui permet d'extraire d'un fichier NASTRAN SOL101 les
%            deplacements sur des points et les efforts sur des elements
%            Les variables en sortie de la fonction sont sauvees dans un
%            binaire Matlab (fichier.mat)
%   cr?? par S.CONIGLIO 17/01/2017
%===========================================================================%
%
% Syntaxe :
% [SUBCASE] = lec_sol101(fichier)
%
% Entree :
%   fichier : chemin du fichier a lire (chaine de caracteres entre ' ')
%
% Sorties :
%  SUBCASE : structure de taille 1*nombre de subcases
%   SUBCASE(i).titre               : titre du subcase i
%   SUBCASE(i).sstitre             : sous-titre du subcase i
%   SUBCASE(i).label               : label du subcase i
%   SUBCASE(i).POINTS(j).num       : numero du point j du subcase i
%   SUBCASE(i).POINTS(j).DISP(k).Mod : deplacement suivant ddl k du point j du subcase i
%   SUBCASE(i).ELEMENTS(j).type    : type de l'element j du subcase i
%   SUBCASE(i).ELEMENTS(j).num     : numero de l'element j du subcase i
%   SUBCASE(i).ELEMENTS(j).ITEM.Mod: effort de l'element j du subcase i
%
% remarque: ne traite que les celas (type 12)
% v3 MV 15/06/07 ajout des elements CELAS11 BAR et ROD
% v4 GT 19/06/07 ajout des elements CONROD
%===========================================================================%
%
% fichier='RB3040_EAS_M11_00_000_BASELINE_Sol101Manvr_01_6sub.pch';
% fichier='RB3040_EAS_M11_00_000_Sol101Manvr_01DP79.pch';
% fichier='RB3040_EAS_M11_00_000_BASELINE_Sol101Manvr_01_Simone.pch';
% fichier='FEM_Level1_08BSC.pch';
% fichier='FEM_Level1_08_M1_000SC.pch';
% fichier='FEM_Level1_08_M3_000SC.pch';
% fichier='FEM_Level2_11BSC.pch';
% fichier='FEM_Level2_11M1SC.pch';
% fichier='FEM_Level2_11M1_000.pch';
% fichier='FEM_Level2_11M3_000.pch';
% fichier='FEM_Level2_11M3SC.pch';
% fichier='FEM_Level1_08.pch';
% fichier='FEM_Level1_08_M1_000.pch';
% fichier='desin_zone_level1-08.dat';
% fichier='complete_design_zone108.dat';
% fichier='level108BASELINE.dat';
% fichier='Level108M1_000.dat';
% fichier='Level2_11baseline.dat';
% fichier='Level108M3_000.dat';
fichier='lin_pert10.dat';
% fichier='FEM_Level2_11.pch';
fid = fopen(fichier,'r');      % ouverture en lecture seule


%
%
compt = 1;
compt_elm = 1;
nbe_point = 0;
nbe_subcase = 0;
nbe_elm = 0;
suite = 1;
count=1;
Bound = 0;
texte   = fgetl(fid);
count = count + 1; LowBound = floor(count/10000)*10000;
if (Bound~=LowBound)
    Bound = LowBound;
    display(['       [' num2str(Bound) ' lines]'])
end

intercoord=[];
%
%
while feof(fid) == 0
    if length(texte)>=length('               O U T P U T   F O R   L O A D   C A S E');
        compt = 1;
        intercount=0;
        if strcmp(texte(1:length('               O U T P U T   F O R   L O A D   C A S E')),'               O U T P U T   F O R   L O A D   C A S E')
            titre= strrep(texte,'               O U T P U T   F O R   L O A D   C A S E','');
            titre=deblank(titre);
            
            for i=1:6
                clear texte
                texte   = fgetl(fid);
                count = count + 1; LowBound = floor(count/10000)*10000;
                if (Bound~=LowBound)
                    Bound = LowBound;
                    display(['       [' num2str(Bound) ' lines]'])
                end
            end
            while ~strcmp(texte(1:length('                                       N O D E   O U T P U T')),'                                       N O D E   O U T P U T')
                for i=1:6
                    clear texte
                    texte   = fgetl(fid);
                    count = count + 1; LowBound = floor(count/10000)*10000;
                    if (Bound~=LowBound)
                        Bound = LowBound;
                        display(['       [' num2str(Bound) ' lines]'])
                    end
                end
            end
            if strcmp(texte(1:length('                                       N O D E   O U T P U T')),'                                       N O D E   O U T P U T')
                nbe_subcase = nbe_subcase+1;
                
                SUBCASE(nbe_subcase).titre = titre;
                SUBCASE(nbe_subcase).sstitre = titre;
                SUBCASE(nbe_subcase).label = titre;
                for i=1:5
                    clear texte
                    texte   = fgetl(fid);
                    count = count + 1; LowBound = floor(count/10000)*10000;
                    if (Bound~=LowBound)
                        Bound = LowBound;
                        display(['       [' num2str(Bound) ' lines]'])
                    end
                end
                while ~isempty(texte)
                    label=strrep(texte,'   THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET ','');
                    for i=1:4
                        clear texte
                        texte   = fgetl(fid);
                        count = count + 1; LowBound = floor(count/10000)*10000;
                        if (Bound~=LowBound)
                            Bound = LowBound;
                            display(['       [' num2str(Bound) ' lines]'])
                        end
                    end
                    texte   = fgetl(fid);
                    count = count + 1; LowBound = floor(count/10000)*10000;
                    if (Bound~=LowBound)
                        Bound = LowBound;
                        display(['       [' num2str(Bound) ' lines]'])
                    end
                    while (isempty(texte)==0)
                        if strcmp(label,'INTERFACE')
                                intercount=intercount+1;
                                intercoord(intercount,:)=[str2double(texte(1:11)),str2double(texte(17+6*12:27+6*12)),str2double(texte(17+7*12:27+7*12)),str2double(texte(17+8*12:27+8*12))];
                            
                        end
                        SUBCASE(nbe_subcase).POINTS(compt).label    = label;
                        SUBCASE(nbe_subcase).POINTS(compt).num    = str2double(texte(1:11));
                        SUBCASE(nbe_subcase).POINTS(compt).DISP(1).Mod = str2double(texte(1+16:11+16));
                        SUBCASE(nbe_subcase).POINTS(compt).DISP(2).Mod = str2double(texte(17+12:27+12));
                        SUBCASE(nbe_subcase).POINTS(compt).DISP(3).Mod = str2double(texte(17+2*12:27+2*12));
                        SUBCASE(nbe_subcase).POINTS(compt).rotation(1).Mod    = str2double(texte(17+3*12:27+3*12));
                        SUBCASE(nbe_subcase).POINTS(compt).rotation(2).Mod    = str2double(texte(17+4*12:27+4*12));
                        SUBCASE(nbe_subcase).POINTS(compt).rotation(3).Mod    = str2double(texte(17+5*12:27+5*12));
                        SUBCASE(nbe_subcase).POINTS(compt).coord(1).Mod    = str2double(texte(17+6*12:27+6*12));
                        SUBCASE(nbe_subcase).POINTS(compt).coord(2).Mod    = str2double(texte(17+7*12:27+7*12));
                        SUBCASE(nbe_subcase).POINTS(compt).coord(3).Mod    = str2double(texte(17+8*12:27+8*12));
                        texte   = fgetl(fid);
                        count = count + 1; LowBound = floor(count/10000)*10000;
                        if (Bound~=LowBound)
                            Bound = LowBound;
                            display(['       [' num2str(Bound) ' lines]'])
                        end
                        compt=compt+1;
                    end
                    for i=1:7
                        clear texte
                        texte   = fgetl(fid);
                        count = count + 1; LowBound = floor(count/10000)*10000;
                        if (Bound~=LowBound)
                            Bound = LowBound;
                            display(['       [' num2str(Bound) ' lines]'])
                        end
                    end
                    if strcmp(texte,' CT: CYLINDRICAL TRANSFORMATION')
                        for i=1:2
                            clear texte
                            texte   = fgetl(fid);
                            count = count + 1; LowBound = floor(count/10000)*10000;
                            if (Bound~=LowBound)
                                Bound = LowBound;
                                display(['       [' num2str(Bound) ' lines]'])
                            end
                        end
                    end
                    for i=1:2
                        clear texte
                        texte   = fgetl(fid);
                        count = count + 1; LowBound = floor(count/10000)*10000;
                        if (Bound~=LowBound)
                            Bound = LowBound;
                            display(['       [' num2str(Bound) ' lines]'])
                        end
                    end
                end
            end
        end
    end
    texte   = fgetl(fid);
    count = count + 1; LowBound = floor(count/10000)*10000;
    if (Bound~=LowBound)
        Bound = LowBound;
        display(['       [' num2str(Bound) ' lines]'])
    end
end
toc
%





%==============================================================================
% Fermeture du fichier
%==============================================================================

fclose(fid);

%==============================================================================
% Sauvegarde dans un fichier binaire
%==============================================================================

pt = findstr(fichier,'.');
eval([ 'save ' fichier '.mat SUBCASE intercoord']);
disp([ 'LE FICHIER ' fichier '.mat A ETE CREE']);
