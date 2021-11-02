%  behavior_master_script_JC
%  play  this for each folder contening "board-DIN-00.dat"
%  "board-DIN-01.dat" ...
%  JCatanese 2017 in  JaegerLab

%%
clear all, close all,
%% To run a loop over folder
MouseID = 'JCVGAT';
AnalysisFolder = 'D:\DATA EMORY\JC_Analysis'; 
% FolderID = dir(['D:\JC_Data\Acute\Awake\' MouseID '*\*day22*\*task*']);
% FolderID = dir(['D:\JC_Analysis\' MouseID '*\*task*']);
FolderID = dir([AnalysisFolder MouseID '*\*taskopto*']);

for nf=1:max(size(FolderID))
    clearvars -EXCEPT FolderID MouseID nf
    cd([FolderID(nf).folder '\' FolderID(nf).name])
    disp([num2str(nf) '/' num2str(max(size(FolderID))) '   MouseID =' MouseID  '   FolderID = ' FolderID(nf).name])
    %% use behavioral events restricted to the task (Start and Stop Manually selected).
%     Get_evt_Restrict_TaskEpoch_JCscript
%     %% Make the behavioral "table" structure containing a binary vec for each basic trial types   
%     Get_TrialsBinVec_TASKrestrict_JCscript
%     %% Count all possible combination of trials types 
    Count_tr_types_JC_Script;
    %% plot1: 4 main trials types   
%     BarPlot_4TrialType_1Sess_BehavOpto_JCscript
%     % plot2 : ALL trials types 
%     BarPlot_AllTrialType_1Sess_BehavOpto_JCscript
 
end
%% Mean over all sessions 
% BarPlot_4TrialType_MeanSTAT_BehavOpto_JCscript
% pub_fig_BehavOpto_BarMeanStat_JCscript
