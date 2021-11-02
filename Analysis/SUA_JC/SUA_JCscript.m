%  SUA_JCscript
% by JC 11-22-2018
%%
clear all, close all,

Mouse = 'JCVGAT';
FolderID = dir(['D:\JC_Analysis\' Mouse '*\*w*d*taskopto*']);
% FolderID = dir(['D:\JC_Analysis\' MouseID '\**\*w*d*taskopto*']);

for nf=1:max(size(FolderID))
    disp([num2str(nf) '/' num2str(max(size(FolderID))) '   MouseID =' Mouse  '   FolderID = ' FolderID(nf).name])
    cd([FolderID(nf).folder '\' FolderID(nf).name]);
    clearvars -except  Mouse FolderID
    close all,

    %% load data
    load('info.mat'); MouseID = info.info_notes.MouseID;  Day=info.info_notes.Day;  sr=info.info_freq_parameters.board_dig_in_sample_rate;
    load('evt.mat');
    load('time.mat');
    load('Ntrial_type.mat');
    
    MouseID
    Day 
    
    delete('*SDF*.jpg')
    delete('*ttest*')
    
    psth_trig_evt = 'Delay'  % center on : 'Delay', % 'APuff' % 'GoCue' % 'Valve' %'Licks'
    psth_trial_type = {'cor'}; % cor = correct ; imp= impulses; Nol = Nolick
    col={'k'};
    SDF_Raster_mlibJCscript
    
    
    psth_trig_evt = 'Delay'  % center on : 'Delay', % 'APuff' % 'GoCue' % 'Valve' %'Licks'
    psth_trial_type = {'cCR', 'cCL'}; % cor = correct ; imp= impulses; Nol = Nolick
    col={['b'] , ['r']};
    SDF_Raster_mlibJCscript
    
end

% STATS TTEST for L-R firing rate
% Stats_plot_LvR_epoch_JC_script



%% list of choices

% SDF_Raster_PSTH_LvRvO_mlibJC_fun('GoCue');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Licks');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'eNO');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'ePL');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'iCR');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'iCL');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'oNO');
% SDF_Raster_PSTH_LvRvO_mlibJC_fun('Delay', 'ocC');










