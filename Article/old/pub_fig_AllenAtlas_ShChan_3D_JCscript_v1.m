% pub_fig_AllenAtlas_ShChan_3D_JCscript 
% save eps figure of the reconstructed electrode placement within the Allen Atlas 
% by JC 2/19/2019

close all; 
load('listcell.mat')
%% Plot 2-D Channels for each single Sessions in color (to indicate which is VM/VL or Others)
NbSess = max(Tcoord.nSess); 
for ii = 1:NbSess
    figure,
    pause(0.2)
    scatter3(Tcoord.AP(Tcoord.nSess==ii), Tcoord.ML(Tcoord.nSess==ii), Tcoord.DV(Tcoord.nSess==ii), 'k')
    hold on, scatter3(Tcoord.AP(Tcoord.VM & Tcoord.nSess==ii), Tcoord.ML(Tcoord.VM & Tcoord.nSess==ii), Tcoord.DV(Tcoord.VM & Tcoord.nSess==ii),'r')
    hold on, scatter3(Tcoord.AP(Tcoord.VL & Tcoord.nSess==ii), Tcoord.ML(Tcoord.VL & Tcoord.nSess==ii), Tcoord.DV(Tcoord.VL & Tcoord.nSess==ii),'b')
    xlabel('AP'), ylabel('ML'), zlabel('DV')
    legend('other', 'VM', 'VL')
    title([ Tcoord.MouseID(min(find(Tcoord.nSess==ii)),:) ' ' Tcoord.Day(min(find(Tcoord.nSess==ii)),:)])
end

%% Plot 3-D ALL Channels ALL sessions in color (to indicate which is VM/VL or Others)
figure, scatter3(Tcoord.AP, Tcoord.ML, Tcoord.DV, 'k')
hold on, scatter3(Tcoord.AP(Tcoord.VM), Tcoord.ML(Tcoord.VM), Tcoord.DV(Tcoord.VM),'r')
hold on, scatter3(Tcoord.AP(Tcoord.VL), Tcoord.ML(Tcoord.VL), Tcoord.DV(Tcoord.VL),'b')
xlabel('AP'), ylabel('ML'), zlabel('DV')
legend('other', 'VM', 'VL')
nVMVL=sum(Tcombo.VMVL);ntot=max(Tcombo.ncell); 
title(['nVM/VL=' num2str(nVMVL) '  ntot=' num2str(ntot)])

%% Plot 2-D Allen Atlas (Coronal or Sagital)
annotation_volume_location = 'C:\Users\JCATANE\Documents\ALLEN BRAIN ATLAS\annotation_volume_10um_by_index.npy';
structure_tree_location = 'C:\Users\JCATANE\Documents\ALLEN BRAIN ATLAS\structure_tree_safe_2017.csv';
template_volume_location = 'C:\Users\JCATANE\Documents\ALLEN BRAIN ATLAS\template_volume_10um.npy';

av = readNPY(annotation_volume_location);
st = loadStructureTree(structure_tree_location);
tv = readNPY(template_volume_location); % tv= template volume = 3D matrice of brain images (10um slices)
size(tv) % ans= 1320         800        1140

for ii = 1:NbSess
    
    Bregma=allenCCFbregma;
    
    AP = -(Tcoord.AP(Tcoord.nSess==ii)*100)+Bregma(1);
    ML = (Tcoord.ML(Tcoord.nSess==ii)*100)+Bregma(3);
    DV = -((Tcoord.DV(Tcoord.nSess==ii)/590)*720*100)+Bregma(2);
    
    if size(unique(AP),1)==1;
        find(strcmp(st.acronym, 'VM'));
        st(646,:)
        figure; im = sliceOutlineWithRegionVec(squeeze(av(unique(AP),:,:)), 646, [1 0 0], gca);
        hold on, sliceOutlineWithRegionVec(squeeze(av(unique(AP),:,:)), 645, [0 0 1], gca);
        hold on; scatter(ML ,DV,10,'g','filled')
        pause(0.2)
    elseif size(unique(ML),1)==1;
        find(strcmp(st.acronym, 'VM'));
        st(646,:)
        figure; im = sliceOutlineWithRegionVec(squeeze(av(:,:,unique(ML))), 646, [1 0 0], gca);
        hold on, sliceOutlineWithRegionVec(squeeze(av(:,:,unique(ML))), 645, [0 0 1], gca);
        hold on; scatter(DV, AP ,15,'g')
        pause(0.2)
    end
    
    
    pause(0.2)
    
Mouse= unique(cellstr(Tcombo.MouseID(Tcombo.nSess==ii,:)));
Day=unique(cellstr(Tcombo.Day(Tcombo.nSess==ii,:)));

    disp(['Saving : D:\JC_Figures\Fig_Article\Histo\Allen_' Mouse{1} '_' Day{1} ])
    saveas(gcf, ['D:\JC_Figures\Fig_Article\Histo\Allen_' Mouse{1} '_' Day{1} ],'eps')
    saveas(gcf, ['D:\JC_Figures\Fig_Article\Histo\Allen_' Mouse{1} '_' Day{1} ],'emf')
    saveas(gcf, ['D:\JC_Figures\Fig_Article\Histo\Allen_' Mouse{1} '_' Day{1} ],'jpg')
    
end