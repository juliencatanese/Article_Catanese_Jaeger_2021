% pub_disp_BrainAreaOpto_JCscript
% Julien Catanese 3-12-2019.

load listcell
%% Display list of Brain area recorded 

annotation_volume_location = 'C:\Users\JCATANE\Documents\ALLEN BRAIN ATLAS\annotation_volume_10um_by_index.npy';
structure_tree_location = 'C:\Users\JCATANE\Documents\ALLEN BRAIN ATLAS\structure_tree_safe_2017.csv';
template_volume_location = 'C:\Users\JCATANE\Documents\ALLEN BRAIN ATLAS\template_volume_10um.npy';

try
    av(1,:);
    st(1,:);
    tv(1,:);
catch
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
    tv = readNPY(template_volume_location); % tv= template volume = 3D matrice of brain images (10um slices)
end

%%
NbArea= max(size(unique(Tcoord.AreaID)));  % 12
AreaUni=  unique(Tcoord.AreaID)
for jj = 1:NbArea
    for iii=1:max(size((st))) ;
        AA= strcmp(st.acronym{iii},AreaUni{jj});
        if AA;
            id{jj}=iii ;
        end;
    end;
    shortname = AreaUni{jj};
    fullname = st.name(id{jj},:);
    NbRecCell = sum(strcmp(Tcoord.AreaID, AreaUni{jj}));
    
    display([ '#'  num2str(NbRecCell) ' ' shortname ' Cells '  fullname{1} ' ( Allen idx = ' num2str(id{jj})  ' ) '])
end

%% display and list of opto Cells area
Ninib  = sum(Tcombo.Opto_inib & (Tcoord.VMVL))
Nexcit = sum(Tcombo.Opto_exct & (Tcoord.VMVL))

disp(['opto inib in VM = ' num2str(Ninib) ])
Tcoord(Tcombo.Opto_inib,:)

disp(['opto Excit in VM = ' num2str(Nexcit) ])
Tcoord(Tcombo.Opto_exct,:)

%% Disp Count Ncell per Sess: VM/VL
VMVL=Tcoord.VMVL
for NN=1:15
    disp(['nVMVLcell= ' num2str(sum(VMVL(Tcoord.nSess==NN))) ' '  Tcoord.MouseID(min(find(Tcoord.nSess==NN)),:) ' '  Tcoord.Day(min(find(Tcoord.nSess==NN)),:)]);
end
%% Disp Count Ncell per Sess: OPTO INHIB
ALL_VMVL_Opto_Inib = sum(Tcombo.VMVL & Tcombo.Opto_inib)
for NN=1:15
    disp(['OPTO_VMVL_inib= ' num2str(sum(Tcombo.nSess==NN & Tcombo.VMVL & Tcombo.Opto_inib)) ' '  Tcombo.MouseID(min(find(Tcombo.nSess==NN)),:) ' '  Tcombo.Day(min(find(Tcombo.nSess==NN)),:)]);
end
%% Disp Count Ncell per Sess: OPTO EXCIT
ALL_VMVL_Opto_Excit = sum(Tcombo.VMVL & Tcombo.Opto_exct)
for NN=1:15
    disp(['OPTO_VMVL_Excit= ' num2str(sum(Tcombo.nSess==NN & Tcombo.VMVL & Tcombo.Opto_exct)) ' '  Tcombo.MouseID(min(find(Tcombo.nSess==NN)),:) ' '  Tcombo.Day(min(find(Tcombo.nSess==NN)),:)]);
end

