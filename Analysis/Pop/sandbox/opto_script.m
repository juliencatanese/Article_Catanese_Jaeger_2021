% Opto_Script
% Written by Julien Catanese 10/27/2018


%% Define Session list
clear all
close all
cd('D:\JC_Analysis');
% SessList = dir(['**/*' '_' 'taskopto' '_' '*']);
% SessList = dir(['**/*taskopto*']);
SessList = dir(['*/*taskopto*']);


NSess= max(size(SessList)) % Number of Sessions
Nol=[]; Nolo=[]; cCR=[]; cCRo=[]; cCL=[]; cCLo=[];

%% loop trhough all Sessions named "taskopto"
for of=1:NSess
    SessID= SessList(of).name
    SessPath=SessList(of).folder;
    cd([SessPath '\' SessID])
    
    %% convert .dat into .mat
    %     dat2mat_JC_Script
    %     info_notes_JC_Script
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
    
    
    %% BEHAVIOR: Table, Count and BarPlot
    
%     Get_Table_trials_behav_JC_Script                    %Make a binnary Table of all trials
%     Count_tr_types_JC_Script                            %Count Separate trial type and save their index
%     BarPlot_AllTrialType_1Sess_BehavOpto_JCscript       %Bar Plot All trial types
%     BarPlot_4TrialType_1Sess_BehavOpto_JCscript         %Bar Plot 4 trials types (ipsi contra ommission impulse)
    
    
    %% EPHYS: Raster, SDF, PSTH
    % SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay');
    % SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'eNO');
    % SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'oNO');
    % SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'ocC');
    
    % SDF_Raster_PSTH_LvRvO_mlibJC_fun('GoCue');
    % SDF_Raster_PSTH_LvRvO_mlibJC_fun('GoCue', 'oNO');
    % SDF_Raster_PSTH_LvRvO_mlibJC_fun('GoCue', 'eNO');
    % SDF_Raster_PSTH_LvRvO_mlibJC_fun('GoCue', 'ocC');
    
    % SDF_Raster_PSTH_LvRvO_mlibJC_fun('Licks');
    % SDF_Raster_PSTH_LvRvO_mlibJC_fun('Licks', 'iCR');
    % SDF_Raster_PSTH_LvRvO_mlibJC_fun('Licks', 'iCL');
    
    % STATS TTEST for L-R firing rate
    
    % Stats_plot_LvR_epoch_JC_script
    
    %%  CLASSIFIER L-R trials based on Firing rate whitin 750ms windows: Compare 3 Epochs
    % linClass_JC_test1
    
     Pop_Seq_DistMatrix_1session_GOcent_corr_mlibJCscript
     ALL_SDFALL = [ALL_SDFALL; sort_sdfall] 
    
end


% figure, plot(NN)

ALL_SDFALL_sorted = sortrows(ALL_SDFALL,'descend');
figure, imagesc(ALL_SDFALL_sorted)
colorbar, 
caxis([0.5 1])
title('ALL SESSIONS')




%% PLOT STAT BEHAVIOR
% BarPlot_4TrialType_MeanSTAT_BehavOpto_JCscript









