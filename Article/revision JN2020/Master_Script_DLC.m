close all, clear all
cd ('D:\DATA EMORY\JC_Analysis')
SAVEFIG_folder =  'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\'
% VGAT14-w14d2
% DATA_folder = 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d2_z4670_VM_taskopto_optopost_G912_180709_vidY_100tr_29cel_10mW'
% MouseID =  'vgat14w14d2'
% load ('Table_DLC_vgat14w14d2.mat')

% VGAT14-w14d8
DATA_folder = 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d8_z4900_VM_taskopto_optopost_G912_180715_vidY_150tr_20cel_13mW_tooDeep'
MouseID =  'vgat14w14d8'
cd (DATA_folder)
load ('vgat14w14d8_DLCresults.mat')

% % VGAT15-w10d8
% DATA_folder = 'D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d7_z4400_VM_taskopto_optopost_G912_180721_vidY_250tr_24cel_11mW'
% MouseIDfull =  'vgat15w10d7'
% cd (DATA_folder)
% load ('Table_DLC_vgat15w10d8.mat')

% VGAT17-w10d5
% DATA_folder = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d5_z4300_VM_taskopto_optopost_EOBD_180728_vidY_200tr_20cel_12mW_Q4'
% MouseID =  'vgat17w10d5'
% cd (DATA_folder)
% load ('vgat17w10d5_DLCresults.mat')

% % VGAT17-w10d7
% DATA_folder = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW'
% MouseID =  'vgat17w10d7'
% cd (DATA_folder)
% load ('vgat17w10d7_DLCresults.mat')


% average over trials (100 frames)
Get_DLC_VAR_Trial_Mat_script

%% OPTIONAL: plot time continuous (raw) DLC-VAR and average (flat)
% plot_continous_DLCvar_script

%% Plot DLC results using accurate trials type selection (corL v corR) or (opto vs omi) ...
plot_save_real_trialtype_average_GoCueCenter

%% TO compare TTL Lick Sensor vs DLC Lick detection
Plot_DLCLick_vs_TTLSensor_align

%% To Analyse single trials
clear all;
close all;
clc;

DATA_folder = 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d8_z4900_VM_taskopto_optopost_G912_180715_vidY_150tr_20cel_13mW_tooDeep'
MouseID =  'vgat14w14d8'

cd (DATA_folder)
load ('vgat14w14d8_DLCresults.mat')
load('Ntrial_type.mat')
load ('evt.mat'); clear evt_delay evt_lick_L evt_lick_R evt_opto evt_puff_L evt_puff_R evt_rwd_L evt_rwd_R;
load ('time.mat');
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);
idx_trial_start = trig_st(1:end-1);
idx_trial_end = idx_trial_start + 4*20000;
% get trigtimes
psth_trig_evt = 'GoCue';
psth_trial_type = {'cor'};
FileLocation= DATA_folder;
[trigtime]=Get_trigtimes(psth_trig_evt, psth_trial_type, FileLocation)
%% Get spks times for 1 trial

FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d8_z4900_VM_taskopto_optopost_G912_180715_vidY_150tr_20cel_13mW_tooDeep'
MouseID =  'vgat14w14d8'
CellID='S3Ch2clu1'
trialtype = 'opt';

Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype)


%%
% for vgat14d8= S3Ch2clu#1

pre=2000
post=2000
NbTtrialType=1;
psth_center_evt=psth_trig_evt;

close all,
CellID='S3Ch2clu1'
load('times_S3Ch2_sub.mat');
CLUST=1;
idx_spk = find(cluster_class(:,1)==CLUST);
spxtimes_ms = cluster_class(idx_spk,2);
spxtimes = spxtimes_ms/10^3;

for ii=1:trial.Nb_correct_L
    figure, hold on,
    
    tr = trial.idx_correct_L(ii)
    idx_tr = [idx_trial_start(tr):1:idx_trial_end(tr)];
    Xt =time(idx_tr); Xt(1);
    Xframe=(Xt-Xt(1))*25;
    
    idx_Spk_tr = find(spxtimes>idx_trial_start(tr)/20000 & spxtimes<idx_trial_end(tr)/20000);
    spxtimes_tr = spxtimes(idx_Spk_tr);
    
    % Compute SDF
    fr=1; %if '1', normalizes to firing rate (Hz); if '0', does nothing (default)
    binsz=1; %ms  bin size of psth (default: 1 ms)
    [psth trialspx] = mpsth(spxtimes_tr, trigtime, 'pre', pre, 'post', post, 'fr',fr, 'binsz', binsz, 'chart',0);
    %PLOT RASTER
    rastmat = zeros(numel(trialspx),pre+1+post);
    timevec = -pre:1:post;
    for i = 1:numel(trialspx);
        rastmat(i,trialspx{i}+pre+1) = 1;
        plot(timevec(rastmat(i,:)~=0),rastmat(i,rastmat(i,:)~=0)*-0.05,'+','Color','r','MarkerSize',10),hold on;
    end
    gca, axis([-pre+10 post+10 0.5 numel(trialspx)+0.5]);
    % PLOT SDF
    ftype = 'Gauss' ;% boxcar, Gauss, exp, exGauss
    w = 40 ;%ms % default was 75ms
    [sdf kernel] = msdf(psth,ftype,w);
    hold on, plot(sdf(:,1), sdf(:,2),'color','k','LineWidth',1 );

    BLtongue= 2*mean(median(allTab.TongueY));
    BLnose= 2*mean(median(allTab.NoseY));
    BLwhisk= 2*mean(median(allTab.WhiskerY));
    
    TongueY_DLC = allTab.TongueY(((tr-1)*100)+1:tr*100);
    hold on, plot(([1:1:100]*40)-1500, (TongueY_DLC/BLtongue)-0.21,'b');
    NoseY_DLC = allTab.NoseY(((tr-1)*100)+1:tr*100);
    hold on, plot(([1:1:100]*40)-1500, (NoseY_DLC/BLnose)-0.21,'g');
    WhiskY_DLC = allTab.WhiskerY(((tr-1)*100)+1:tr*100);
    hold on, plot(([1:1:100]*40)-1500, (WhiskY_DLC/BLwhisk)-0.21,'y');
    
    legend('Spks', 'PSTH', 'TongueY', 'NoseY', 'WhiskerY'  )
    xlabel(['time from ' psth_center_evt ' start ( ms)'], 'FontSize',11,'FontWeight', 'normal');
    ylim([-0.1 0.2+max(sdf(:,2))])
    xlim([-750 750])
    title(['trial ' num2str(tr)], 'FontSize',11,'FontWeight', 'normal');
    hold on, line([-750 -750], [-2 2], 'Color','k','LineStyle','--','LineWidth',1 );
    hold on, line([0 0], [-2 2], 'Color','k','LineStyle','--','LineWidth',1 );
    hold on, plot(sdf(:,1),ones(1,max(size(sdf)))*-0.05,'r');
    title([MouseID ' ' CellID ' trial ' num2str(tr)], 'FontSize',11,'FontWeight', 'normal');
    legend('Spks', 'PSTH', 'TongueY', 'NoseY', 'WhiskerY'  )
    saveas(gcf, ['Fig_' MouseID '_SingleTrial#' num2str(tr) '_' CellID ],'png')
    
end
