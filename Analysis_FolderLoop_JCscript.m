% Analysis_FolderLoop_JCscript
%  written by Julien Catanese in 2017 in JaegerLab
% last updated: 12/07/2018

%%
close all,
cd 'D:\JC_Analysis\'
%% To run a loop over folder
MouseID = 'JCVGAT17';
FolderID = dir(['D:\JC_Analysis\' MouseID '*\*w10d8*task*']);

for nf=1:max(size(FolderID))
    disp([num2str(nf) '/' num2str(max(size(FolderID))) '   MouseID =' MouseID  '   FolderID = ' FolderID(nf).name])
    cd([FolderID(nf).folder '\' FolderID(nf).name])
    
    %% convert .dat into .mat
    %     dat2mat_JC_Script
    %     info_notes_JC_Script
    %     remap_native2custom_JC_Script
    %     %% delete the .dat
    %     delete *.dat
    %     delete *A-0*_raw.mat
    
    %     %% Automatic Spike Sorting using Wave_clus
    %     Artifact_removal_Refmean_JC_Script
    %     % Artifact_removal_1ChanRef_JC_fun(Ref_Chan)
    %     Auto_wave_clus_JCscript % threshold spike detection (5std) + clustering process (see Quiroga paper).
    %     disp('CHECK MANUAL CLUSTERING')
    delete('*SDF*')
    close all
    SDF_Raster_PSTH_LvRvO_mlibJC_fun('GoCue');
    close all
    dashedline=0;
    psth_center_evtID = 'Licks'
    SDF_Raster_PSTH_Correct_Impulse_Opto_mlib_JCscript
    close all
    dashedline=1;
    psth_center_evtID = 'GoCue'
    SDF_Raster_PSTH_Correct_Impulse_Opto_mlib_JCscript
    close all
    
end
dd
%%
close all
violet=[0.5 0.2 0.7]
Sh=2
Ch=7

% % SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, psth_trig_evt , psth_trial_type, col, pre, post)
CLUST=1;
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'GoCue' , {'cCL','cCR'}, {'r','b'}, 1600, 1600)
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'GoCue' , {'cor','imp','eNO'}, {'b',[0.5 0.2 0.7],'c'}, 1600, 1600)
% SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'Licks' , {'cor','imp'}, {'b',[0.5 0.2 0.7],'c'}, 1600, 1600)

%
CLUST=2;
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'GoCue' , {'cCL','cCR'}, {'r','b'}, 1600, 1600)
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'GoCue' , {'cor','imp','eNO'}, {'b',[0.5 0.2 0.7],'c'}, 1600, 1600)
% SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'Licks' , {'cor','imp'}, {'b',[0.5 0.2 0.7],'c'},1600, 1600)
%

CLUST=3;
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'GoCue' , {'cCL','cCR'}, {'r','b'}, 1600, 1600)
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'GoCue' , {'cor','imp','eNO'}, {'b',[0.5 0.2 0.7],'c'}, 1600, 1600)
% SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'Licks' , {'cor','imp'}, {'b',[0.5 0.2 0.7],'c'}, 1600, 1600)
%%

CLUST=4;
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'GoCue' , {'cCL','cCR'}, {'r','b'}, 1600, 1600)
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'GoCue' , {'cor','imp','eNO'}, {'b',[0.5 0.2 0.7],'c'}, 1600, 1600)
% SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'Licks' , {'cor','imp'}, {'b',[0.5 0.2 0.7],'c'}, 1600, 1600)

%%
CLUST=5;
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'GoCue' , {'cCL','cCR'}, {'r','b'}, 1600, 1600)
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'GoCue' , {'cor','imp','eNO'}, {'b',[0.5 0.2 0.7],'c'}, 1600, 1600)
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'Licks' , {'cor','imp'}, {'b',[0.5 0.2 0.7],'c'}, 1600, 1600)

close all
%%

CLUST=6;
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'GoCue' , {'cCL','cCR'}, {'r','b'}, 1600, 1600)
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'GoCue' , {'cor','imp','eNO'}, {'b',[0.5 0.2 0.7],'c'}, 1600, 1600)
SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, 'Licks' , {'cor','imp'}, {'b',[0.5 0.2 0.7],'c'}, 1600, 1600)

%

