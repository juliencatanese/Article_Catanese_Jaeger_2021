close all, clear all, clc,
cd ('D:\DATA EMORY\JC_Analysis')

%% Video Analysis
Plot_DLC_AverageAcrossMice
%% Single trial Analysis
clc, close all
% Mouse#3 Cell#1
FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d8_z4900_VM_taskopto_optopost_G912_180715_vidY_150tr_20cel_13mW_tooDeep';
MouseID =  'vgat14w14d8';
CellID='S3Ch2clu1';
cd(FileLocation);
% trtype='cor'; Plot_DLC_results_JCScript;
% PLOT_RASTER_PSTH_JCscript
trialtype = 'cor';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);
trialtype = 'imp';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);
trialtype = 'omi';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);


%%
% Mouse#4 Cell#2
FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d7_z4400_VM_taskopto_optopost_G912_180721_vidY_250tr_24cel_11mW';
MouseID =  'vgat15w10d7';
CellID='S1Ch1clu4';
cd(FileLocation);
% trtype='cor'; Plot_DLC_results_JCScript;
% PLOT_RASTER_PSTH_JCscript
trialtype = 'cor';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);
trialtype = 'imp';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);
trialtype = 'omi';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);

%% Mouse#5 Cell#3
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW';
MouseID =  'vgat17w10d7';
CellID='S4Ch5clu1';
cd(FileLocation);
% trtype='cor'; Plot_DLC_results_JCScript;
% PLOT_RASTER_PSTH_JCscript
trialtype = 'cor';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);
trialtype = 'imp';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);
trialtype = 'omi';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);

%% Example for Figure3 IMP OMI
% cell#74  (vgat12w11d5S3Ch5clu#02)
% cell#198 (vgat14w14d8S3Ch2clu#01)
% cell#377 (vgat17w10d4S2Ch3clu#01)

%% Cell#74 = vgat12w11d5S3Ch5clu#02
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT12\vgat12_w11d5_z4300_VM_taskopto_optopost_CCAE_180609_vidY_100tr_42cel_10mW_3otr';
MouseID =  'vgat12w11d5';
CellID='S3Ch5clu2'; Sh=str2num(CellID(2));Ch=str2num(CellID(5)); CLUST=str2num(CellID(9));
psth_trial_type={'cor'}; col={'k'};pre=1000; post=1000; K=1;  GaussSmooth=40;
psth_center_evt='Licks';
pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'png')
saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'emf')

%% cell#377 = vgat17w10d4S2Ch3clu#01
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d4_z4200_VM_taskopto_optopost_CCAE_180726_vidY_075tr_41cel_10mW_bl_4otr';
MouseID =  'vgat17w10d4';
CellID='S2Ch3clu1'; Sh=str2num(CellID(2));Ch=str2num(CellID(5)); CLUST=str2num(CellID(9));
psth_trial_type={'cor'}; col={'k'};pre=1000; post=1000; K=1;  GaussSmooth=40;
psth_center_evt='Licks';
pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'png')
saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'emf')

%% ADDITIONAL EXAMPLES
%% Mouse#4 Cell#4
FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d7_z4400_VM_taskopto_optopost_G912_180721_vidY_250tr_24cel_11mW';
MouseID =  'vgat15w10d7';
CellID='S2Ch8clu1';
cd(FileLocation);
% trtype='cor'; Plot_DLC_results_JCScript;
% PLOT_RASTER_PSTH_JCscript
trialtype = 'cor';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);
trialtype = 'imp';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);
trialtype = 'omi';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);

%% Mouse#5 Cell#5
FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d5_z4300_VM_taskopto_optopost_EOBD_180728_vidY_200tr_20cel_12mW_Q4';
MouseID =  'vgat17w10d5';
CellID='S1Ch1clu2';
cd(FileLocation);
% trtype='cor'; Plot_DLC_results_JCScript;
% PLOT_RASTER_PSTH_JCscript
trialtype = 'cor';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);
trialtype = 'imp';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);
trialtype = 'omi';
Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype);




%% MULTI-CELL SINGLE TRIALS FIGURE (7cells in 5cor, 2imp, 2omi trials )
% load Tfig1_VMopto.mat
% Tfig1_VMopto(Tfig1_VMopto.nSess==14 & Tfig1_VMopto.VMVL==1 & Tfig1_VMopto.Opto_inib ==1 ,:)
clear all, clc,
MouseID =  'vgat17w10d7';
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW';
cd(FileLocation);
trialtype = 'cor';
Xlimite = [-750 750]; 
close all ;
%% 7 cells from Mouse#5 ALL VM/VAL and opto inhib and Ramping
CellID='S3Ch3clu1'; [sdf1, Xsdf, tridx] = Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype, Xlimite); close all
CellID='S2Ch6clu1'; [sdf2, Xsdf, tridx] = Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype, Xlimite);close all
CellID='S2Ch7clu1'; [sdf3, Xsdf, tridx] = Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype, Xlimite);close all
CellID='S3Ch5clu1'; [sdf4, Xsdf, tridx] = Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype, Xlimite);close all
CellID='S4Ch5clu1'; [sdf5, Xsdf, tridx] = Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype, Xlimite);close all
CellID='S4Ch6clu2'; [sdf6, Xsdf, tridx] = Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype, Xlimite);close all
CellID='S4Ch7clu2'; [sdf7, Xsdf, tridx] = Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype, Xlimite);close all
% SAVING
cd ('D:\DATA EMORY\JC_Analysis')
save(['7cells_sdf_alltr_' trialtype '_fig5_5STrials'],'sdf1','sdf2','sdf3','sdf4','sdf5','sdf6','sdf7', 'tridx', 'Xsdf')

%%
clear all, cd ('D:\DATA EMORY\JC_Analysis')

trialtype = 'omi'
load(['7cells_sdf_alltr_' trialtype '_fig5_5STrials'],'sdf1','sdf2','sdf3','sdf4','sdf5','sdf6','sdf7', 'tridx', 'Xsdf')

if trialtype == 'cor';
    TrID= [4; 65; 117; 147; 196]
elseif trialtype == 'imp';
    TrID= [15; 145]
elseif trialtype == 'omi';
    TrID= [44; 91]
end

for it= 1:max(size(TrID))
    idxTr(it) = find(tridx==TrID(it))
end

sdf1 = sdf1(:,idxTr)';
sdf2 = sdf2(:,idxTr)';
sdf3 = sdf3(:,idxTr)';
sdf4 = sdf4(:,idxTr)';
sdf5 = sdf5(:,idxTr)';
sdf6 = sdf6(:,idxTr)';
sdf7 = sdf7(:,idxTr)';

SDF7cell=[]; SDF7cell_MEAN=[];
for it= 1:max(size(TrID))
    SDF7cell{it} =              [sdf1(it,:); sdf2(it,:); sdf3(it,:); sdf4(it,:); sdf5(it,:); sdf6(it,:); sdf7(it,:)]
    SDF7cell_MEAN(it,:) = nanmean( [sdf1(it,:); sdf2(it,:); sdf3(it,:); sdf4(it,:); sdf5(it,:); sdf6(it,:); sdf7(it,:)])
end


% %% plot all cells all trials
% close all, clc,
% figure, hold on, plot(Xsdf, (SDF7cell{1})), legend('1','2','3','4','5','6','7', 'Location','best')
% figure, hold on, plot(Xsdf, (SDF7cell{2})), legend('1','2','3','4','5','6','7', 'Location','best')
% % figure, hold on, plot(Xsdf, (SDF7cell{3})), legend('1','2','3','4','5','6','7', 'Location','best')
% % figure, hold on, plot(Xsdf, (SDF7cell{4})), legend('1','2','3','4','5','6','7', 'Location','best')
% % figure, hold on, plot(Xsdf, (SDF7cell{5})), legend('1','2','3','4','5','6','7', 'Location','best')

%% plot cell group average for each trials
% close all
clc
figure, hold on,
Xlimite= [-750 750]
subplot(6,1,1), plot(Xsdf, nanmean(SDF7cell{1})), xlim([Xlimite])
subplot(6,1,2), plot(Xsdf, nanmean(SDF7cell{2})), xlim([Xlimite])
% subplot(6,1,3), plot(Xsdf, nanmean(SDF7cell{3})), xlim([Xlimite])
% subplot(6,1,4), plot(Xsdf, nanmean(SDF7cell{4})), xlim([Xlimite])
% subplot(6,1,5), plot(Xsdf, nanmean(SDF7cell{5})), xlim([Xlimite])
subplot(6,1,6), plot(Xsdf, nanmean(SDF7cell_MEAN), 'r', 'LineWidth',2);
xlim([Xlimite])
% legend('1','2','3','4','5', 'Location','best')
title([trialtype ' trials average (7cells)'] )
% hold on, line([-750 -750], [-1 1], 'Color','k','LineStyle','--','LineWidth',1 );
% hold on, line([0 0], [-1 1], 'Color','k','LineStyle','--','LineWidth',1 );
% ylim([-0.01 0.2])
% hold on, plot(sdf(:,1),ones(1,max(size(sdf)))*-0.05,'r');
%%

FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW';
MouseID =  'vgat17w10d7';
cd(FileLocation);
trtype='cor';

CellID='S1Ch6clu1'; PLOT_RASTER_PSTH_JCscript; sdf_mean1= sdf_mean; sdf_sem1= sdf_sem; 
CellID='S1Ch8clu2'; PLOT_RASTER_PSTH_JCscript; sdf_mean2= sdf_mean; sdf_sem1= sdf_sem; 
CellID='S2Ch3clu2'; PLOT_RASTER_PSTH_JCscript; sdf_mean3= sdf_mean; sdf_sem1= sdf_sem; 
CellID='S2Ch5clu1'; PLOT_RASTER_PSTH_JCscript; sdf_mean4= sdf_mean; sdf_sem1= sdf_sem; 
CellID='S2Ch6clu1'; PLOT_RASTER_PSTH_JCscript; sdf_mean5= sdf_mean; sdf_sem1= sdf_sem; 
CellID='S2Ch7clu1'; PLOT_RASTER_PSTH_JCscript; sdf_mean6= sdf_mean; sdf_sem1= sdf_sem; 
CellID='S3Ch3clu1'; PLOT_RASTER_PSTH_JCscript; sdf_mean6= sdf_mean; sdf_sem1= sdf_sem; 
CellID='S3Ch5clu1'; PLOT_RASTER_PSTH_JCscript; sdf_mean7= sdf_mean; sdf_sem1= sdf_sem; 
CellID='S4Ch5clu1'; PLOT_RASTER_PSTH_JCscript; sdf_mean8= sdf_mean; sdf_sem1= sdf_sem; 
CellID='S4Ch6clu2'; PLOT_RASTER_PSTH_JCscript; sdf_mean9= sdf_mean; sdf_sem1= sdf_sem; 
CellID='S4Ch7clu2'; PLOT_RASTER_PSTH_JCscript; sdf_mean10= sdf_mean; sdf_sem= sdf_sem; 

%% 7 cells PSTH grabd Average over all cor trial
%% A- With 40ms smooth 
close all, 
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW';
MouseID =  'vgat17w10d7';
cd(FileLocation);
trtype='cor';

CellID='S2Ch6clu1'; PLOT_RASTER_PSTH_JCscript; SMA40(1,:)= sdf_mean; SMA40_sem(1,:)= sdf_sem; 
CellID='S2Ch7clu1'; PLOT_RASTER_PSTH_JCscript; SMA40(2,:)= sdf_mean; SMA40_sem(2,:)= sdf_sem; 
CellID='S3Ch3clu1'; PLOT_RASTER_PSTH_JCscript; SMA40(3,:)= sdf_mean; SMA40_sem(3,:)= sdf_sem; 
CellID='S3Ch5clu1'; PLOT_RASTER_PSTH_JCscript; SMA40(4,:)= sdf_mean; SMA40_sem(4,:)= sdf_sem; 
CellID='S4Ch5clu1'; PLOT_RASTER_PSTH_JCscript; SMA40(5,:)= sdf_mean; SMA40_sem(5,:)= sdf_sem; 
CellID='S4Ch6clu2'; PLOT_RASTER_PSTH_JCscript; SMA40(6,:)= sdf_mean; SMA40_sem(6,:)= sdf_sem; 
CellID='S4Ch7clu2'; PLOT_RASTER_PSTH_JCscript; SMA40(7,:)= sdf_mean; SMA40_sem(7,:)= sdf_sem; 
%
SDF_7mean = mean(SMA40);
SDF_7sem = mean(SMA40_sem);

%
figure, 
plot([-pre:1:post],SDF_7mean, 'r')
hold on, 
plotshaded([-pre:1:post],[SDF_7mean-SDF_7sem; SDF_7mean+SDF_7sem], 'r')
xlim(Xlimite)
%% 7 cells PSTH grabd Average over all cor trial
%% B- With 75ms smooth (SMA)
close all

load('SMA_cor_GoCue545.mat')
pre=parfigUsed.pre;
post=parfigUsed.post;
IDXcellID = [474 476 479 480 489 492 494]
IDXcellID = [489]
IDXcellID = [198]

AA = SMA(IDXcellID,:)
BB = SSemA(IDXcellID,:)
AAA = mean(AA)
BBB = mean(BB)
SSS = std(AA,1)/sqrt(7)


Xgocue = [-pre:1:post]
figure, hold on 
plot(Xgocue, AA'), 
title('7 cells separated')
xlim(Xlimite)

figure, hold on
plot(Xgocue, AAA)
title('Average 7 cells ')
xlim(Xlimite)

figure, hold on
plot(Xgocue, AAA,'r','LineWidth',2)
plotshaded(Xgocue, [AAA-SSS; AAA+SSS],'r')
title('Average 7 cells ')
xlim(Xlimite)

figure, hold on
plot(Xgocue, AAA,'r','LineWidth',2)
plotshaded(Xgocue, [AAA-2*BBB; AAA+2*BBB],'r')
title('Average 7 cells ')
xlim(Xlimite)

%% CLASSIFIER FIGURE
% parameters
close all;
mypath = 'D:\DATA EMORY\JC_Analysis'
cd(mypath)
parfig.pathinit = mypath
SaveFigFolder = 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig7_classifier'
clearvars -except mypath parfig SaveFigFolder
parfig.ControlShuffle=0;
rng('shuffle');
load('listcell3.mat'); 
load('Tfig1_VMopto.mat');

parfig.Nrepeat = 100;
parfig.typeClass='FoldXVal' ; %     typeClass='%HoldOut'
parfig.learner= 'logistic';%, 'svm'}
parfig.Nfold = 10;

parfig.center_evt = 'GoCue';% center ('Delay' ; 'GoCue' ;'APuff'; 'Licks')
parfig.epoch = 'delay'

if parfig.epoch  =='delay'
    parfig.pre = 750 ; % define how much time before zero (in ms)
    parfig.post= 0;
elseif parfig.epoch == 'Apuff'
    parfig.pre = 1500 ; % define how much time before zero (in ms)
    parfig.post= -750;
end

%% Compute All raw traces per sessions (very long to recompute)
% 11 session with 15 cells NOL
% 12 sessions with 20 cells VMVL
celltype = 'VMVL'
idxcell = Tfig1_VMopto.VMVL
parfig.minNbcell = 20;
parfig.ControlShuffle=0;

trial_type_1 = 'imp'
trial_type_2 = 'omi'
pub_fig7_Classifier_JCscript % very long



%% Figure6: PLOT Average for distinct population
% PRE-LICK epoch, Lick Centered, cor vs imp
% PLOTING Mean +/- 1SEM (shaded)
close all 
SaveFigFolder = 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig7_classifier'
clearvars -except mypath parfig Tfig1_VMopto listcell SaveFigFolder jj

parfig.center_evt = 'GOcue';    % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.epoch = 'delay'
% 
% figure, hold on,
for ii=1:3
    if ii==1
        celltype = 'VMVL'; LegendName{ii} = celltype;
        parfig.minNbcell = 20;
        col = 'k'
    elseif ii==2
        celltype = 'ShuffleControl'; LegendName{ii} = celltype;
        parfig.minNbcell = 20;
        col = 'r'
    elseif ii==3
        celltype = 'SvTh-'; LegendName{ii} = celltype;
        parfig.minNbcell = 9;
        col = 'c'
    end
    
    for jj=1:3
        if jj==1
            f1= figure(1), hold on,
            trial_type_1 = 'cor'
            trial_type_2 = 'imp';
        elseif jj==2
            f2=figure(2), hold on,
            trial_type_1 = 'cor'
            trial_type_2 = 'omi';
        elseif jj==3
            f3=figure(3), hold on,
            trial_type_1 = 'imp'
            trial_type_2 = 'omi';
        end
        
        
        load(['Class10Fold_Mall_' trial_type_1 'v' trial_type_2  '_' parfig.epoch '_' parfig.center_evt '_' celltype '_minCell' num2str(parfig.minNbcell) '.mat'], 'Mall')
        
        K = 1;
        Mm=mean(Mall)
        NSessFinal = size(Mall,1)
        Mstd = std(Mall)/sqrt(NSessFinal)
        XX=[1:1:parfig.minNbcell];
        plot( XX , Mm, col, 'LineWidth', 3)
        plotshaded(XX, [Mm-K*Mstd ; Mm+K*Mstd], col)
        title(['Mean Class10Fold ' trial_type_1 ' vs ' trial_type_2  ' ' parfig.epoch  ' ' parfig.center_evt ' '  celltype ])


    end
end

hold on,
plot( XX , ones(1,parfig.minNbcell)/2, 'k--')
ylim([0.2 0.8]);
title(['Mean Class10Fold ' trial_type_1 ' vs ' trial_type_2  ' ' parfig.epoch  ' ' parfig.center_evt ' '  celltype ])
legend(LegendName{1}, '', LegendName{2}', '',LegendName{3},'')

%% FIGURE OPTO RASTER
close all, 
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW';
cd(FileLocation);
MouseID =  'vgat17w10d7';
ListCellSameSess = {'S1Ch6clu1'; 'S1Ch8clu2';...
    'S2Ch3clu2'; 'S2Ch5clu1'; 'S2Ch6clu1';...
    'S2Ch7clu1'; 'S3Ch3clu1'; 'S3Ch5clu1';...
    'S4Ch5clu1'; 'S4Ch6clu2'; 'S4Ch7clu2'}

for ii=1:max(size(ListCellSameSess))
CellID=ListCellSameSess{ii} 

Sh=str2num(CellID(2));Ch=str2num(CellID(5)); CLUST=str2num(CellID(9));
psth_center_evt='GoCue'; 
psth_trial_type={'cor','opt'}; col={'k','c'};
pre=1200; post=1200; K=1;  
GaussSmooth=40; 
pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)

end
