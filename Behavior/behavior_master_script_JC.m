%  behavior_master_script_JC
%  play  this for each folder contening "board-DIN-00.dat"
%  "board-DIN-01.dat" ...
%  JCatanese 2017 in  JaegerLab

%%
clear all, close all,
%% To run a loop over folder
MouseID = 'JCVGAT';
% FolderID = dir(['D:\JC_Data\Acute\Awake\' MouseID '*\*day22*\*task*']);
FolderID = dir(['D:\JC_Analysis\' MouseID '*\*task*']);

for nf=1:max(size(FolderID))
    clearvars -EXCEPT FolderID MouseID nf
    cd([FolderID(nf).folder '\' FolderID(nf).name])
    disp([num2str(nf) '/' num2str(max(size(FolderID))) '   MouseID =' MouseID  '   FolderID = ' FolderID(nf).name])
   
    %% Count and BarPlot Separate trial type and save idx
    Get_Table_trials_behav_JC_Script;
    %%
    Count_tr_types_JC_Script;
    %%
    BarPlot_4TrialType_1Sess_BehavOpto_JCscript
    %%
    BarPlot_AllTrialType_1Sess_BehavOpto_JCscript
end

BarPlot_4TrialType_MeanSTAT_BehavOpto_JCscrip
