% pub_fig_AllenAtlas_ShChan_3D_JCscript
% save eps figure of the reconstructed electrode placement within the Allen Atlas
% by JC 2/19/2019

% clear all
FigSave_ON = 0; 
myfolder = 'C:\Users\catan\Documents\EMORY\JC_Analysis\AllenBrainAtlas\'

annotation_volume_location = [myfolder 'annotation_volume_10um_by_index.npy'];
structure_tree_location = [myfolder  'structure_tree_safe_2017.csv'];
template_volume_location = [myfolder  'template_volume_10um.npy'];

load ('Tfig1_VMopto.mat');  


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

close all
nSess=max(Tfig1_VMopto.nSess);
for ii = 15%1:1% nSess
    
    Bregma=allenCCFbregma;
    
    AP = -(Tfig1_VMopto.AP(Tfig1_VMopto.nSess==ii)*100)+Bregma(1);
    ML = (Tfig1_VMopto.ML(Tfig1_VMopto.nSess==ii)*100)+Bregma(3);
    DV = -((Tfig1_VMopto.DV(Tfig1_VMopto.nSess==ii)/590)*720*100)+Bregma(2);
    
    APoi = -(Tfig1_VMopto.AP(Tfig1_VMopto.nSess==ii & Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_inib)*100)+Bregma(1);
    MLoi = (Tfig1_VMopto.ML(Tfig1_VMopto.nSess==ii & Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_inib)*100)+Bregma(3);
    DVoi = -((Tfig1_VMopto.DV(Tfig1_VMopto.nSess==ii & Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_inib)/590)*720*100)+Bregma(2);
    
    APoe = -(Tfig1_VMopto.AP(Tfig1_VMopto.nSess==ii & Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_exct)*100)+Bregma(1);
    MLoe = (Tfig1_VMopto.ML(Tfig1_VMopto.nSess==ii & Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_exct)*100)+Bregma(3);
    DVoe = -((Tfig1_VMopto.DV(Tfig1_VMopto.nSess==ii & Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_exct)/590)*720*100)+Bregma(2);
    
    APvm = -(Tfig1_VMopto.AP(Tfig1_VMopto.nSess==ii & Tfig1_VMopto.VMVL)*100)+Bregma(1);
    MLvm = (Tfig1_VMopto.ML(Tfig1_VMopto.nSess==ii & Tfig1_VMopto.VMVL)*100)+Bregma(3);
    DVvm = -((Tfig1_VMopto.DV(Tfig1_VMopto.nSess==ii & Tfig1_VMopto.VMVL)/590)*720*100)+Bregma(2);
    
    APvl = -(Tfig1_VMopto.AP(Tfig1_VMopto.nSess==ii & Tfig1_VMopto.VMVL)*100)+Bregma(1);
    MLvl = (Tfig1_VMopto.ML(Tfig1_VMopto.nSess==ii & Tfig1_VMopto.VMVL)*100)+Bregma(3);
    DVvl = -((Tfig1_VMopto.DV(Tfig1_VMopto.nSess==ii & Tfig1_VMopto.VMVL)/590)*720*100)+Bregma(2);
    
    
    MLall= []; 
    DVall= []; 
    step=(0.1/590)*720*100+Bregma(2); 
    for ns=1:8
        DVall= [DVall round(max(DV)-(step*(ns-1)))];
        MLall= [MLall unique(ML)];
    end

    DVall = [DVall; DVall; DVall; DVall]
    
    Mouse= unique(cellstr(Tfig1_VMopto.MouseID(Tfig1_VMopto.nSess==ii,:)));
    Day=unique(cellstr(Tfig1_VMopto.Day(Tfig1_VMopto.nSess==ii,:)));
    idxVM= find(strcmp(st.acronym, 'VM')); % idxVM = 646
    idxVL= find(strcmp(st.acronym, 'VAL')); % idxVL = 645
    st(idxVM,:)
    
    
    APU=unique(AP)
    NbAPU = size(APU,1)
    if   NbAPU==1
        
        figure,  sliceOutlineWithRegionVec(squeeze(av(unique(AP),:,:)), idxVM, 'r' , gca);
        hold on, sliceOutlineWithRegionVec(squeeze(av(unique(AP),:,:)), idxVL, 'b', gca);
        for ii=1:1 
            hold on, scatter(MLall(ii,:), DVall(ii,:),12,'k','filled')
%             hold on, scatter(MLall(ii,:), DVall(ii,:),12,'k')
        end
        hold on, scatter(MLvl ,DVvl, 6,'w', 'filled'), hold on, scatter(MLvl ,DVvl, 12,'k')
        hold on, scatter(MLvm ,DVvm, 6,'w', 'filled'), , hold on, scatter(MLvm ,DVvm, 12,'k')
        hold on, scatter(MLoi ,DVoi, 12,'c') 
        hold on, scatter(MLoe ,DVoe, 12,'m') 
        title([  Mouse{1} ' ' Day{1}  '(AP=' num2str(round(unique(AP)-Bregma(1))/100) ')' ])
        pause(0.2)
        
        
    else
        figure, 
        for nm=NbAPU:-1:1
            idxAPU = find(AP==APU(nm))
            subplot(2,2, nm)
            im = sliceOutlineWithRegionVec(squeeze(av(round(APU(nm)),:,:)), 646, [1 0 0], gca);
            hold on, sliceOutlineWithRegionVec(squeeze(av(round(APU(nm)),:,:)), 645, [0 0 1], gca);
            hold on; scatter(ML(idxAPU) ,DV(idxAPU),8,'g')
            title([  Mouse{1} ' ' Day{1} '(AP=' num2str(round(APU(nm)-Bregma(1))/100) ')' ])
        end
        
        pause(0.2)
    end
    
    if FigSave_ON==1
        pause(0.2)
        disp(['Saving : D:\JC_Figures\Fig_Article\AllenAtlas\Allen_' Mouse{1} '_' Day{1} ])
        saveas(gcf, ['D:\JC_Figures\Fig_Article\AllenAtlas\Allen_' Mouse{1} '_' Day{1} ],'eps')
        saveas(gcf, ['D:\JC_Figures\Fig_Article\AllenAtlas\Allen_' Mouse{1} '_' Day{1} ],'emf')
        saveas(gcf, ['D:\JC_Figures\Fig_Article\AllenAtlas\Allen_' Mouse{1} '_' Day{1} ],'jpg')
        pause(0.2)
    end
end

%% 3D PLOT With VM and VL separated

for nc= 1:Tfig1_VMopto.ncell(end) %545
  VAL(nc,1)=  logical(strcmp(Tfig1_VMopto.AreaID{nc},'VAL'));   
end

X0= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL));
Y0= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL));
Z0= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL));

X1= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & VAL));
Y1= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & VAL));
Z1= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & VAL));

X2= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_inib));
Y2= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_inib));
Z2= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_inib));

X3= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_exct));
Y3= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_exct));
Z3= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_exct));

figure, 
scatter3(X0(:),Y0(:),Z0(:),'k','filled'), view(-60,60)
hold on 
scatter3(X1(:),Y1(:),Z1(:),'k','filled'), view(-60,60)
hold on, 
scatter3(X2(:),Y2(:),Z2(:),'c','filled'), view(-60,60)
hold on, 
scatter3(X3(:),Y3(:),Z3(:),'m','filled'), view(-60,60)
xlabel('AP')
ylabel('ML')
zlabel('DV')
legend ('VM','VAL','Th-','Th+')

