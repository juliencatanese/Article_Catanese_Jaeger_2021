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
trtype='cor'
Plot_DLC_results_JCfun(FileLocation, MouseID, trtype);
%%
% trialtype = 'cor';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'imp';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'omi';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'opO';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'opI';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'opC';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
cd(FileLocation); 
% close all;
% Raster SDF Shaded
% GoCue 
% Sh=str2num(CellID(2));Ch=str2num(CellID(5)); CLUST=str2num(CellID(9));
% psth_trial_type={'cor'}; col={'k'};pre=1000; post=1000; K=1;  GaussSmooth=40; 
% psth_center_evt='GoCue'; 
% pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'png')
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'emf')
% % Licks 
% psth_center_evt='Licks'; 
% pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'png')
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'emf')
% Plot Average DLC results 
Get_Abs_Average_3DLCvar_JCScript
plot_Abs_AverageLR_3DLCvar_JCScript
plot_Abs_Average_3DLCvar_JCScript2
%%
% Mouse#4 Cell#2
FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d7_z4400_VM_taskopto_optopost_G912_180721_vidY_250tr_24cel_11mW';
MouseID =  'vgat15w10d7';
CellID='S1Ch1clu4';
Plot_DLC_results_JCfun(FileLocation, MouseID);
% trialtype = 'cor';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'imp';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'omi';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'opO';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'opI';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'opC';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
%% 
cd(FileLocation);
% close all;
% Raster SDF Shaded
% GoCue 
% Sh=str2num(CellID(2));Ch=str2num(CellID(5)); CLUST=str2num(CellID(9)); 
% psth_trial_type={'cor'}; col={'k'};pre=1000; post=1000; K=1;  GaussSmooth=40; 
% psth_center_evt='GoCue'; 
% pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'png')
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'emf')
% Licks 
% psth_center_evt='Licks'; 
% pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'png')
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'emf')
% Plot Average DLC results 
% Get_Abs_Average_3DLCvar_JCScript
% plot_Abs_AverageLR_3DLCvar_JCScript
% plot_Abs_Average_3DLCvar_JCScript2


%% Mouse#5 Cell#3
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW';
MouseID =  'vgat17w10d7';
CellID='S4Ch5clu1';
cd(FileLocation); 

Plot_DLC_results_JCfun(FileLocation, MouseID);
% trialtype = 'cor';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'imp';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'omi';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'opO';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'opI';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'opC';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
%% 
% close all;
% Raster SDF Shaded
%GoCue 
% Sh=str2num(CellID(2));Ch=str2num(CellID(5)); CLUST=str2num(CellID(9)); 
% psth_trial_type={'cor'}; col={'k'};pre=1000; post=1000; K=1;  GaussSmooth=40; 
% psth_center_evt='GoCue'; 
% pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'png')
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'emf')
% Licks 
% psth_center_evt='Licks'; 
% pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'png')
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'emf')
% Plot Average DLC results 
Get_Abs_Average_3DLCvar_JCScript
plot_Abs_AverageLR_3DLCvar_JCScript
plot_Abs_Average_3DLCvar_JCScript2

dddd


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
% Mouse#4 Cell#4
FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d7_z4400_VM_taskopto_optopost_G912_180721_vidY_250tr_24cel_11mW';
MouseID =  'vgat15w10d7';
% Plot_DLC_results_JCfun(FileLocation, MouseID);
CellID='S2Ch8clu1';
% Raster SDF Shaded
%GoCue 
cd(FileLocation); close all;
Sh=str2num(CellID(2));Ch=str2num(CellID(5)); CLUST=str2num(CellID(9)); 
psth_trial_type={'cor'}; col={'k'};pre=1000; post=1000; K=1;  GaussSmooth=40; 
psth_center_evt='GoCue'; 
pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'png')
saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'emf')
% Licks 
psth_center_evt='Licks'; 
pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'png')
saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'emf')
% 
% trialtype = 'cor';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'imp';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'omi';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'opO';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'opI';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% trialtype = 'opC';
% Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
% % Plot Average DLC results 
%%


%%
% Mouse#5 Cell#5
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW';
MouseID =  'vgat17w10d7';
CellID='S4Ch6clu2';
trialtype = 'cor';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)
trialtype = 'imp';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)
trialtype = 'omi';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)
trialtype = 'opO';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)
trialtype = 'opI';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)
trialtype = 'opC';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)


% Mouse#5 Cell#6
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW';
MouseID =  'vgat17w10d7';
CellID='S3Ch3clu1';
trialtype = 'cor';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
trialtype = 'imp';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
trialtype = 'omi';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
trialtype = 'opO';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
trialtype = 'opI';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);
trialtype = 'opC';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype);

% Mouse#5 Cell#7
FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d5_z4300_VM_taskopto_optopost_EOBD_180728_vidY_200tr_20cel_12mW_Q4';
MouseID =  'vgat17w10d5';
% Plot_DLC_results_JCfun(FileLocation, MouseID);
CellID='S1Ch1clu2';
trialtype = 'cor';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)
trialtype = 'imp';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)
trialtype = 'omi';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)
trialtype = 'opO';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)
trialtype = 'opI';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)
trialtype = 'opC';
Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)


