function [SDF_all, XSDF, tridx] = Plot_SingleTrial_raster_JCfun(MouseID, CellID, trialtype, Xlimite)
% function Plot_SingleTrial_raster_JCfun(FileLocation, MouseID, CellID, trialtype, Xlimite)
% Julien Catanese 10/21/2020
%%
close all,
% parameters
psth_center_evt = 'GoCue';
pre=2000 %ms before center_evt
post=2000 %ms after center evt
psth_trial_type = {trialtype};
NbTtrialType=1;
CLUST= str2num(CellID(end));

% load data and events
% cd (FileLocation); mkdir('.\SingleTrial')
load ([ MouseID '_DLCresults.mat']); % video data (DLC, 25Hz)
load('Ntrial_type.mat'); % trial data (TTL, 20KHz)
load ('evt.mat', 'evt_trial'); % events data (TTL, 20KHz)
load ('time.mat'); % times data (20KHz)
load(['times_' CellID(1:5) '_sub.mat']); %spks data (20KHz)

%% Get trigtimes
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);
idx_trial_start = trig_st(1:end-1);
idx_trial_end = idx_trial_start + 4*20000;
FileLocation = pwd;
[trigtime]=Get_trigtimes(psth_center_evt, psth_trial_type, FileLocation)

% Get spks times for 1 trial
idx_spk = find(cluster_class(:,1)==CLUST);
spxtimes_ms = cluster_class(idx_spk,2);
spxtimes = spxtimes_ms/10^3;

% Get trials to analyses
if trialtype == 'cor'
    Ntr= trial.Nb_correct_L + trial.Nb_correct_R
    tridx = sort([trial.idx_correct_L trial.idx_correct_R])
elseif trialtype == 'omi'
    Ntr= trial.Nb_NoLick
    tridx = sort([trial.idx_NoLick])
elseif trialtype == 'imp'
    Ntr= trial.Nb_errorDelay_PL_CL + trial.Nb_errorDelay_PL_CR + trial.Nb_errorDelay_PR_CL + trial.Nb_errorDelay_PR_CR
    tridx = sort([trial.idx_errorDelay_PL_CL trial.idx_errorDelay_PR_CR trial.idx_errorDelay_PR_CL trial.idx_errorDelay_PL_CR])
elseif trialtype == 'opC'
    Ntr= trial.Nb_correct_L_opto + trial.Nb_correct_R_opto
    tridx = sort([trial.idx_correct_L_opto trial.idx_correct_R_opto])
elseif trialtype == 'opO'
    Ntr= trial.Nb_NoLick_opto
    tridx = sort([trial.idx_NoLick_opto])
elseif trialtype == 'opI'
    Ntr= trial.Nb_errorDelay_PL_CL_opto + trial.Nb_errorDelay_PL_CR_opto + trial.Nb_errorDelay_PR_CL_opto + trial.Nb_errorDelay_PR_CR_opto
    tridx = sort([trial.idx_errorDelay_PL_CL_opto trial.idx_errorDelay_PR_CR_opto trial.idx_errorDelay_PR_CL_opto trial.idx_errorDelay_PL_CR_opto])
end

% loop across trials
SDF_all=[]; 
for ii=1:Ntr
    figure, hold on,
    
    tr =  tridx(ii)
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
%     size(sdf)
%     sdf
    SDF_all = [SDF_all sdf(:,2)/(max(sdf(750:2250,2)))];
    XSDF= sdf(:,1); 
    
    hold on, plot(sdf(:,1), sdf(:,2)/(5*max(sdf(:,2))),'color','k','LineWidth',1.5 );
    
    BLtongue= 2*mean(median(allTab.TongueY));
    BLnose= 2*mean(median(allTab.NoseY));
    BLwhisk= 2*mean(median(allTab.WhiskerY));
    
    NoseY_DLC = allTab.NoseY(((tr-1)*100)+1:tr*100);
    hold on, plot(([1:1:100]*40)-1500, (NoseY_DLC/BLnose)-0.20,'c');

    WhiskX_DLC = allTab.WhiskerX(((tr-1)*100)+1:tr*100);
    hold on, plot(([1:1:100]*40)-1500, (WhiskX_DLC/BLwhisk),'m');
    
    TongueY_DLC = allTab.TongueY(((tr-1)*100)+1:tr*100);
    hold on, plot(([1:1:100]*40)-1500, (TongueY_DLC/BLtongue)-0.30,'g','LineWidth',1.5);

    hold on, line([-750 -750], [-2 2], 'Color','k','LineStyle','--','LineWidth',1 );
    hold on, line([0 0], [-2 2], 'Color','k','LineStyle','--','LineWidth',1 );
    hold on, plot(sdf(:,1),ones(1,max(size(sdf)))*-0.05,'r');
        
    title([trialtype ' trial ' num2str(tr) ' ' MouseID ' ' CellID ' trial ' num2str(tr)], 'FontSize',11,'FontWeight', 'normal');
    xlabel(['time from ' psth_center_evt ' start ( ms)'], 'FontSize',11,'FontWeight', 'normal');
    ylim([-0.1 0.50])
    xlim([Xlimite])
    set(gca,'YTick',[])
    legend('Spikes', 'PSTH',  'NoseY', 'WhiskerX','TongueY', 'Location','northwest')
    saveas(gcf, ['.\SingleTrial\Fig_'  MouseID  '_SingleTrial#' num2str(tr) '_' CellID '_' trialtype],'png')
    mkdir(['.\SingleTrial\' trialtype ])
    saveas(gcf, ['.\SingleTrial\' trialtype '\Fig_'  MouseID '_' trialtype '_SingleTrial#' num2str(tr) '_' CellID ],'png')
    saveas(gcf, ['.\SingleTrial\'  trialtype '\Fig_' MouseID '_' trialtype '_SingleTrial#' num2str(tr) '_' CellID  ],'emf')
        
end