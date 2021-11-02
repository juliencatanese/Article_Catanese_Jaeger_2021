% pub_table2_Tcoord_JCscript
% create and save a table ('T2') that contain all cells informations
% Fields: ncell, nMouse, nDay, nChan, nShank, AP, ML, DV, VM, VL.
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'T2')
% written by Julien Catanese 12/05/2018
% last updated JC 2/7/2019.
% last updated JC 4/12/2019. Add Tcombo. reduce size Tcoord. 
% 
clearvars -except mypath parfig
load('listcell.mat');

AP_measured=[]; ML_measured=[]; DV_measured=[];
AP_corr=[]; ML_corr=[]; DV_corr=[]; C=[]; B=[]; A=[];
T1=[]; T1= listcell; T1=removevars(T1,[5,6,7,8,9])

%% Add nSess (Session number from 1:nSess instead of the string Day)
listSess=unique(listcell.SessID)
nSess = max(size(listSess))
for ii= 1:nSess
    A(:,ii)= strcmp(listcell.SessID, listSess{ii});
end
B=A.*[1:1:nSess];
C=sum(B')'; % vector with Session Number from 1 to 11.
T1= addvars(T1, C, 'Before', 5, 'NewVariableNames', {'nSess'});
T1(1,:);

%% Add nMouse (Mice number from 1:nMouse instead of the string MouseID)
ncell=T1.ncell(end);
C=[]; B=[]; A=[];
for ii=1:ncell;
    MID{ii}= T1.MouseID(ii,:);
end
listMice=unique(MID');
nMice = max(size(listMice))
for ii= 1:nMice
    A(:,ii)= strcmp(MID, listMice{ii});
end
B=A.*[1:1:nMice];
C=sum(B')'; % vector with Session Number from 1 to 11.
T1= addvars(T1, C, 'Before', 6, 'NewVariableNames', {'nMouse'});
T1(1,:)

%% Manual Coordinates (based on Stereotaxic Measurement and Histo reconstruction and Ephys Markers)
% 'Sh1M_' = The coordinate in the table corespond to channel #8 (always tip) on Shank #1 (the most Medial M_)


SessionNAME = ['11_w10d4'; '12_w11d5'; '14_w14d2'; '14_w14d5'; '14_w14d7'; '14_w14d8'; '15_w08d5'; '15_w10d3'; '15_w10d7'; '15_w10d8'; '17_w10d3'; '17_w10d4'; '17_w10d5'; '17_w10d7'; '17_w10d8'  ];
Ch_Origin   = ['S1Ch8M_' ; 'S2Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S4Ch8MA' ; 'S1Ch8LA' ;  'S1Ch8LP'; 'S1Ch8M_' ; 'S2Ch8_L' ; 'S1Ch4_Q' ; 'S1Ch8M_' ;  'S1Ch8M_'  ];
Ch_Orient   = ['ML'      ;  'AP'     ;  'AP'     ;  'AP'     ;  'AP'     ;  'AP'     ;  'AP'     ; 'DG'      ;  'DG'     ;  'DG'     ;  'ML'     ;  'ML'     ;  'Q4'     ;    'ML'   ;   'ML'      ];
ElectrodeID = ['CC4F_n4S';'CCAE_o3S' ;'G912_o4S' ;'G912_o4S' ;'CC4F_n4S' ;'G912_o4S' ;'G912_o4S' ; 'CAED_n4S';'G912_o4S' ;'G912_o4S' ;'CAED_n4S' ; 'CCAE_o3S'; 'EOBD_o1S'; 'G912_o4S'; 'G912_o4S'  ];
AP_measured = [1.35      ;  1.58     ;  1.10     ;  1.55     ;  1.50     ;  1.55     ;  1.90     ;  1.20     ;   1.2     ;   1.5     ;  1.00     ;  1.00     ;   1.3     ;     1.45  ;     1.2     ];
ML_measured = [0.35      ;  0.90     ;  0.75     ;  0.90     ;  0.80     ;  0.85     ;  0.70     ;  0.95     ;   0.8     ;  0.75     ;  0.5      ;  1.10     ;   0.65    ;     0.65  ;     0.65    ];
DV_measured = [4.35      ;  4.30     ;  4.65     ;  4.30     ;  4.35     ;  4.90     ;  4.70     ;  4.25     ;   4.4     ;   4.6     ;  4.30     ;  4.20     ;   4.30    ;     4.30  ;     4.5     ];
APAdj       = [+0.40     ;  +0.00    ;  +0.15    ;  +0.15    ;  +0.15    ;  +0.15    ;  -0.20    ;  -0.20    ;  -0.20    ;  -0.20    ;  +0.20    ;  +0.00    ;   -0.15   ;     +0.0  ;     +0.0    ];    % due to Bending of the Headpost or electrode.
MLAdj       = [+0.10     ;  +0.00    ;  +0.00    ;  +0.00    ;  +0.00    ;  +0.00    ;  +0.20    ;  +0.15    ;  +0.10    ;  +0.25    ;  +0.15    ;  +0.15    ;   +0.15   ;     +0.15 ;    +0.15    ];    % due to Bending of the Headpost or electrode.
DVAdj       = [+0.20     ;  +0.20    ;  +0.10    ;  +0.25    ;  +0.25    ;  +0.25    ;  +0.00    ;  +0.05    ;  +0.05    ;  +0.15    ;  +0.10    ;  +0.20    ;   +0.20   ;     +0.20 ;    +0.020   ]; % due to increasing cortical lesion after several implants and some departure from the bones.

DVAdj=DVAdj-0.20; % due to difference between touch cortex and touch water. 


AP_corr = AP_measured + APAdj;
ML_corr = ML_measured + MLAdj;
DV_corr = DV_measured + DVAdj ;

Tst= table( Ch_Origin, Ch_Orient, AP_corr, ML_corr, DV_corr, listSess)

%%
nShank = str2num(T1.ChanID(:,2));
nChan  = str2num(T1.ChanID(:,5));
T0=[]; T0=T1;
T2=addvars(T0,nShank,nChan); T2(1,:)

%%
AP=[]; ML=[];DV=[];
for ii=1:ncell;
    %   coln =;
    Orig=Tst.Ch_Origin(T1.nSess(ii),:);
    ShOrig=str2num(Orig(2));
    ChOrig=str2num(Orig(5));
    Orient=Orig(6:7);
    if ShOrig~=4
        if Orient=='M_'
            AP(ii) = Tst.AP_corr(T1.nSess(ii));
            ML(ii) = Tst.ML_corr(T1.nSess(ii))+(0.2*(-ShOrig+T2.nShank(ii)));
            DV(ii)  = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
            
        elseif Orient=='_L'
            AP(ii) = Tst.AP_corr(T1.nSess(ii));
            ML(ii) = Tst.ML_corr(T1.nSess(ii))-(0.2*(-ShOrig+T2.nShank(ii)));
            DV(ii)  = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
            
        elseif Orient=='_P'
            AP(ii) = Tst.AP_corr(T1.nSess(ii))-(0.2*(-ShOrig+T2.nShank(ii)));
            ML(ii) = Tst.ML_corr(T1.nSess(ii));
            DV(ii)  = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
            
        elseif Orient=='_Q'
            AP(ii) = Tst.AP_corr(T1.nSess(ii));
            ML(ii) = Tst.ML_corr(T1.nSess(ii));
            DV(ii) = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
            
        elseif Orient=='LP'
            AP(ii) = Tst.AP_corr(T1.nSess(ii))-(((600*cos(pi/6))/3000)*(-ShOrig+T2.nShank(ii)));
            ML(ii) = Tst.ML_corr(T1.nSess(ii))-(((600*sin(pi/6))/3000)*(-ShOrig+T2.nShank(ii)));
            DV(ii) = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
        
        elseif Orient=='LA'
            AP(ii) = Tst.AP_corr(T1.nSess(ii))+(((600*cos(pi/6))/3000)*(-ShOrig+T2.nShank(ii)));
            ML(ii) = Tst.ML_corr(T1.nSess(ii))-(((600*sin(pi/6))/3000)*(-ShOrig+T2.nShank(ii)));
            DV(ii) = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
        
        else
            lvkjsdlvjsdlvkj
        end
        
    elseif ShOrig==4
        if Orient=='MA'
            AP(ii) = Tst.AP_corr(T1.nSess(ii))+(((600*cos(pi/6))/3000)*(ShOrig-T2.nShank(ii)));
            ML(ii) = Tst.ML_corr(T1.nSess(ii))-(((600*sin(pi/6))/3000)*(ShOrig-T2.nShank(ii)));
            DV(ii)  = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
            
        end
    end
end
AP=-AP'; ML=-ML'; DV=-DV';
%% Add Coordinates to Table
nShank = str2num(T1.ChanID(:,2));
nChan  = str2num(T1.ChanID(:,5));
T0=[]; T0=T1; T2=[];
T2=addvars(T0, nShank, nChan, AP, ML, DV); T2(1,:)


%%
% pub_fig_AllenAtlas_ShChan_3D_JCscript2
%% Plot 2-D Allen Atlas (Coronal or Sagital)
% myfolder= 'C:\Users\JCATANE\Documents\ALLEN BRAIN ATLAS\'
% myfolder = 'C:\Users\Julien\Documents\WORK\AllenBrainAtlas\'
myfolder = 'C:\Users\catan\Documents\EMORY\JC_Analysis\AllenBrainAtlas\'

annotation_volume_location = [myfolder 'annotation_volume_10um_by_index.npy'];
structure_tree_location = [myfolder  'structure_tree_safe_2017.csv'];
template_volume_location = [myfolder  'template_volume_10um.npy'];

try
    av(1,:);
    st(1,:);
    tv(1,:);
catch
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
    tv = readNPY(template_volume_location); % tv= template volume = 3D matrice of brain images (10um slices)
end
% size(tv) % ans= 1320         800        1140
%% plot 
close all
for ii = 1:5%nSess
    
    Bregma=allenCCFbregma;
    
    AP = -(T2.AP(T2.nSess==ii)*100)+Bregma(1);
    ML = (T2.ML(T2.nSess==ii)*100)+Bregma(3);
    DV = -((T2.DV(T2.nSess==ii)/590)*720*100)+Bregma(2);
    
    Mouse= unique(cellstr(T1.MouseID(T1.nSess==ii,:)));
    Day=unique(cellstr(T1.Day(T1.nSess==ii,:)));
    
    if size(unique(AP),1)==1;
        find(strcmp(st.acronym, 'VM'));
        st(646,:)
        figure; im = sliceOutlineWithRegionVec(squeeze(av(unique(AP),:,:)), 646, [1 0 0], gca);
        hold on, sliceOutlineWithRegionVec(squeeze(av(unique(AP),:,:)), 645, [0 0 1], gca);
        hold on; scatter(ML ,DV,8,'g')
        title([  Mouse{1} ' ' Day{1}  '(AP=' num2str(round(unique(AP)-Bregma(1))/100) ')' ])
        pause(0.2)
    elseif size(unique(ML),1)==1;
        find(strcmp(st.acronym, 'VM'));
        st(646,:)
        figure; im = sliceOutlineWithRegionVec(squeeze(av(:,:,unique(ML))), 646, [1 0 0], gca);
        hold on, sliceOutlineWithRegionVec(squeeze(av(:,:,unique(ML))), 645, [0 0 1], gca);
        hold on; scatter(DV, AP ,5,'g')
        title([  Mouse{1} ' ' Day{1} '(ML=' num2str(round(unique(ML)-Bregma(3))/100) ')' ])
        pause(0.2)
    else
        find(strcmp(st.acronym, 'VM'));
        figure;
        APU=unique(AP)
        NbAPU = size(APU,1)
        
        for nm=1:NbAPU
            idxAPU = find(AP==APU(nm))
            subplot(2,2, nm)
            im = sliceOutlineWithRegionVec(squeeze(av(round(APU(nm)),:,:)), 646, [1 0 0], gca);
            hold on, sliceOutlineWithRegionVec(squeeze(av(round(APU(nm)),:,:)), 645, [0 0 1], gca);
            hold on; scatter(ML(idxAPU) ,DV(idxAPU),8,'g')
            title([  Mouse{1} ' ' Day{1} '(AP=' num2str(round(APU(nm)-Bregma(1))/100) ')' ])
        end
        
        pause(0.2)
    end
    
    if parfig.FigSave_ON==1
        pause(0.2)
        disp(['Saving : D:\JC_Figures\Fig_Article\Histo\Allen_' Mouse{1} '_' Day{1} ])
        saveas(gcf, ['D:\JC_Figures\Fig_Article\Histo\Allen_' Mouse{1} '_' Day{1} ],'fig')
        saveas(gcf, ['D:\JC_Figures\Fig_Article\Histo\Allen_' Mouse{1} '_' Day{1} ],'emf')
        saveas(gcf, ['D:\JC_Figures\Fig_Article\Histo\Allen_' Mouse{1} '_' Day{1} ],'jpg')
        pause(0.2)
    end
end
%%

VM=[]; VL=[]; AreaIDall=[];ZI = []; VPM=[];
for nc=1:ncell
    idxName = av(-round((T2.AP(nc)*100))+Bregma(1), -round((T2.DV(nc)/590)*720*100)+Bregma(2), round((T2.ML(nc)*100)+Bregma(3)));
    AreaID = st.acronym(idxName) ;
    AreaIDall{nc} = AreaID{1};
    VM = [VM; strcmp(AreaID{1},'VM')];
    VL = [VL; strcmp(AreaID{1},'VAL')];
    ZI = [ZI; strcmp(AreaID{1},'ZI')];
    VPM = [VPM; strcmp(AreaID{1},'VPM')];
    disp([ 'Cell #' num2str(nc)  ' is ' num2str(AreaID{1}) ])
end
%% display number of VMVL cells
nbVM= sum(VM)
nbVL= sum(VL)
nbVMVL=nbVM+nbVL

%% Save Tcoord
AreaID=[]; AreaID=AreaIDall;
AreaID = AreaID';
VMVL = VM + VL;
ncell = T1.ncell;
nSess = T1.nSess;
nMouse = T1.nMouse;

Tcoord = table(ncell, nMouse, nSess, VMVL, AreaID) 
disp('Tcoord done but no saved')

Tfig1_VMopto = addvars(T2, VMVL, AreaID)
save([mypath '\listcell.mat'],'Tfig1_VMopto', '-append')
disp('Tfig1_VMopto saved')

Tcombo=[]; 
Tcombo=listcell(:,1:5); AreaID=Tcoord.AreaID; VMVL=Tcoord.VMVL; 
ncell =  Tcoord.ncell; nSess =  Tcoord.nSess; nMouse = Tcoord.nMouse;
Tcombo = addvars(Tcombo, nMouse, nSess ,VMVL, AreaID);
save([mypath '\Tcombo.mat'],'Tcombo'); Tcombo(1,:)
disp('Tcombo SAVED')


