% clear all 

annotation_volume_location = 'C:\Users\JCATANE\Documents\ALLEN BRAIN ATLAS\annotation_volume_10um_by_index.npy';
structure_tree_location = 'C:\Users\JCATANE\Documents\ALLEN BRAIN ATLAS\structure_tree_safe_2017.csv';
template_volume_location = 'C:\Users\JCATANE\Documents\ALLEN BRAIN ATLAS\template_volume_10um.npy';

av = readNPY(annotation_volume_location);
st = loadStructureTree(structure_tree_location);
tv = readNPY(template_volume_location); % tv= template volume = 3D matrice of brain images (10um slices) 
size(tv) % ans= 1320         800        1140
%%
close all 

find(strcmp(st.acronym, 'VM'))
st(646,:)
coord_Ch1 = [1.2 (3.6/590)*720 -1.2]*100+allenCCFbregma % coord = [AP=660   DV=360   ML=450] 
coord_Ch8 = [1.2 (4.2/590)*720 -1.2]*100+allenCCFbregma % coord = [AP=660   DV=420   ML=450] 

figure; im = sliceOutlineWithRegionVec(squeeze(av(660,:,:)), 646, [0 1 0], gca);
hold on; im = sliceOutlineWithRegionVec(squeeze(av(660,:,:)), 645, [0 0 1], gca);

% figure; im = sliceOutlineWithRegionVec(squeeze(av(:,:,450)), 646, [1 0 0], gca);
hold on; plot(430,linspace(coord_Ch1(2),coord_Ch8(2),8),'r','Marker','.','MarkerSize',10,'MarkerFaceColor','r')
hold on; plot(450,linspace(coord_Ch1(2),coord_Ch8(2),8),'r','Marker','.','MarkerSize',10,'MarkerFaceColor','r')
hold on; plot(470,linspace(coord_Ch1(2),coord_Ch8(2),8),'r','Marker','.','MarkerSize',10,'MarkerFaceColor','r')
hold on; plot(490,linspace(coord_Ch1(2),coord_Ch8(2),8),'r','Marker','.','MarkerSize',10,'MarkerFaceColor','r')
% hold on; plot(350,linspace(coord_Ch1(2),coord_Ch8(2),8),'Marker','.','MarkerSize',10,'MarkerFaceColor','r')
% hold on; plot(550,linspace(coord_Ch1(2),coord_Ch8(2),8),'Marker','.','MarkerSize',10,'MarkerFaceColor','r')
% hold on; plot(450,linspace(coord_Ch1(2),coord_Ch8(2),8),'Marker','.','MarkerSize',10,'MarkerFaceColor','r')

% hold on; plot(coord_Ch1,'bo')
% hold on; plot(450, 0:10:80, 'bo')
%%
% figure; imagesc(squeeze(tv(660,:,:)))
% axis off
% colormap gray
% axis equal
% hold on; plot(450,360,'ro')
% size(tv)