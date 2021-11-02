% JC_AnalysisDJlab_Master_Script
%  behavior_master_script_JC
%  play  this for each folder contening "board-DIN-00.dat"
%  "board-DIN-01.dat" ...
%  written by Julien Catanese in 2017 in JaegerLab
% last updated: 10/08/2018

%%
clear all, close all,
%% To run a loop over folder
MouseID = 'JCVGAT11';
% FolderID = dir(['D:\JC_Data\Acute\Awake\' MouseID '*\*day22*\*task*']);
FolderID = dir(['D:\JC_Analysis\' MouseID '*\*w10d4*VM_task*']);
% FolderID = dir(['D:\JC_Data\Acute\Awake\JCVGAT12_ChR2_Acute-VM-S4Ch32o\w11d3_rec\z4300_task_180607_183846'])
for nf=1:max(size(FolderID))
    disp([num2str(nf) '/' num2str(max(size(FolderID))) '   MouseID =' MouseID  '   FolderID = ' FolderID(nf).name])
    cd([FolderID(nf).folder '\' FolderID(nf).name])
    
    %% convert .dat into .mat
    dat2mat_JC_Script
    info_notes_JC_Script 
%     remap_native2custom_JC_Script
%     %% delete the .dat
%     delete *.dat
%     delete *A-0*_raw.mat
%     
%     %% Automatic Spike Sorting using Wave_clus
%     Artifact_removal_Refmean_JC_Script 
%     % Artifact_removal_1ChanRef_JC_fun(Ref_Chan)
%     Auto_wave_clus_JCscript % threshold spike detection (5std) + clustering process (see Quiroga paper). 
    disp('CHECK MANUAL CLUSTERING')
   
    %% plot the behavioral raster trial/trials with licks    
%         plot_task1_LickRASTER_JC_Script
    
    % Count and BarPlot Separate trial type and save idx 
%         Get_Table_trials_behav_JC_Script
%         Count_tr_types_JC_Script
%         Plot_bar_trial_type_JC_Script
%     %% NEED TO DO MANUAL CLUSTERING

    
    %Plot RASTER
%     plot_RASTER_Spk_test4_JC_Script

    %Plot SDF_Raster_PSTH center on Delay
    
    % SDF_Raster_PSTH_LvRvO_mlibJC_fun('APuff');
% % function [psth trialspx] = SDF_Raster_PSTH_LvRvO_mlibJC_fun(psth_center_evtID, Other_trial_type)
% % INPUT1 : evt_center = 
% % 'Delay'; 'APuff'; 'GoCue'; 'Licks'; 'Valve'
% % INPUT2 : Other_trial_type =    
% % 'iCR' = impulse Choice Right (or Left: 'iCL')
% % 'oCR' = Opto correct Choice Right (or Left: 'oCL')
% % 'ocC' = Opto correct Choices (L+R)
% % 'oNO' = Opto Not Licked. 

% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'eNO');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'ePL');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'iCR');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'oNO');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'ocC');
% 
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Licks');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Licks', 'iCR');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Licks', 'iCL');


% SDF_Raster_PSTH_LvRvO_mlibJC_fun('GoCue');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('GoCue', 'oNO');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('GoCue', 'eNO');
% % SDF_Raster_PSTH_LvRvO_mlibJC_fun('GoCue', 'ocC');
% 
% SDF_Raster_PSTH_Correct_Impulse_Opto_mlib_JCscript
% STATS TTEST for L-R firing rate 
% 
% Stats_plot_LvR_epoch_JC_script
% 
%  CLASSIFIER L-R trials based on Firing rate whitin 750ms windows: Compare 3 Epochs
% linClass_JC_test1
% Classifier_logistic_10Xval_cor_imp_JCscript
% Classifier_logistic_10Xval_Left_Right_JCscript


end






