function [sdf_mean, sdf_sem, nbtrial]= pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_trig_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
% function pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_trig_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
% Default INPUT Paramnters
%     psth_trig_evt =  'Delay'; % psth center on Delay start; can also use 'GoCue'
%     psth_trial_type = {'cor'}; % All correct trials; to spearate Left Right: {'cCL', 'cCR'} for Left Right
%     col={'k'} ; or col={'r', 'b'}
%     pre=1500;
%     post=1500;
%  by JC 11/22/2018 last updated 10/12/2020

%% Define NbTtrialType
try
    K
    pre
    post
    col
    psth_trig_evt
    psth_trial_type
    NbTtrialType = max(size(psth_trial_type))
    disp('Success input')
catch
    K=1
    disp('ATTENTION: USING DEFAULT PARAM')
    psth_trig_evt =  'Delay' % center for SDF
    psth_trial_type = {'cor'}
    NbTtrialType = max(size(psth_trial_type));
    col={'k'}
    pre=1500
    post=1500
    disp('ATTENTION: USING DEFAULT PARAM')
end

%% load
load('info.mat'); MouseID = info.info_notes.MouseID;  Day=info.info_notes.Day;  sr=info.info_freq_parameters.board_dig_in_sample_rate;
load('evt.mat');
load('time.mat');
load('Ntrial_type.mat');

%% GET trigtimes : times of the psth_trig_evt is sec   (1 vector of times)
% define trial idx
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);

idx_trial_start = trig_end(1:end-1); %-(1*sr);
idx_trial_end = trig_st(2:end) ;     %-(1*sr);

Ntrials= trial.Ntrial;

% define common events
evt_valve = evt_rwd_L + evt_rwd_R;
evt_lick = evt_lick_L + evt_lick_R;
evt_puff =evt_puff_L + evt_puff_R;
evt_prelick = evt_delay + evt_puff;

%% define trigtimes for each trials
trigtimes=[];
for tr=1:Ntrials
    time_tr = time(idx_trial_start(tr):idx_trial_end(tr));
    if psth_trig_evt=='APuff'
        puff_tr=  evt_puff(idx_trial_start(tr):idx_trial_end(tr));
        idx_puff_st = find(diff(puff_tr)>0); % start of the puff
        time_puff_st = time_tr(idx_puff_st);
        trigtimes = [trigtimes time_puff_st]; % in sec
        
    elseif psth_trig_evt=='Delay'
        delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
        idx_delay_st = find(diff(delay_tr)>0);   % start of the delay
        time_delay_st = time_tr(idx_delay_st);
        trigtimes = [trigtimes time_delay_st]; % in sec
        
    elseif psth_trig_evt=='GoCue'
        delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
        idx_GO_st = find(diff(delay_tr)<0); % end of the delay
        time_GO_st = time_tr(idx_GO_st);
        trigtimes = [trigtimes time_GO_st]; % in sec
        
    elseif psth_trig_evt=='Licks'
        lick_tr=  evt_lick(idx_trial_start(tr):idx_trial_end(tr));
        idx_lick_st = min(find(diff(lick_tr)>0));% start of the first lick
        if isempty(idx_lick_st) % replace by Lick evt by GOcue evt for centering during NOlick trial
            delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
            idx_GO_st = find(diff(delay_tr)<0)
            idx_lick_st=idx_GO_st;
        end
        time_lick_st = time_tr(idx_lick_st);
        trigtimes = [trigtimes time_lick_st]; % in sec
        
    else
        disp(' WRONG INPUT ARG pick one: psth_trig_evt = GoCue, APuff, Delay, Licks or Valve' )
    end
end

%%
trigtimes_all=trigtimes';
trigtimes_cCL = trigtimes(trial.idx_correct_L);
trigtimes_cCR = trigtimes(trial.idx_correct_R);
trigtimes_iCL = trigtimes([trial.idx_errorDelay_PL_CL trial.idx_errorDelay_PR_CL]);
trigtimes_iCR = trigtimes([trial.idx_errorDelay_PL_CR trial.idx_errorDelay_PR_CR]);
trigtimes_oCR = trigtimes([trial.idx_correct_R_opto]);
trigtimes_oCL = trigtimes([trial.idx_correct_L_opto]);
trigtimes_ocC = trigtimes([trial.idx_correct_L_opto trial.idx_correct_R_opto]);
trigtimes_oNO = trigtimes([trial.idx_NoLick_opto]);
trigtimes_eNO = trigtimes([trial.idx_NoLick]);
trigtimes_ePL = trigtimes([trial.idx_errorResp_PL]);
trigtimes_ePR = trigtimes([trial.idx_errorResp_PR]);

trigtimes_cor = sort([trigtimes_cCR trigtimes_cCL]);
trigtimes_imp = sort([trigtimes_iCL trigtimes_iCR]);
trigtimes_opt = sort([trigtimes_oNO trigtimes_ocC]);


%% GET spxtimes = SPIKE times in Sec (1 vector of times)
load(['times_S' num2str(Sh) 'Ch' num2str(Ch) '_sub.mat']); %  e.g.: load('times_S2Ch6_man.mat')
ChanID = ['S' num2str(Sh) 'Ch' num2str(Ch) ]

idx_spk = find(cluster_class(:,1)==CLUST);
spxtimes_ms = cluster_class(idx_spk,2);
spxtimes = spxtimes_ms/10^3;        % convert time to sec

figure,
NbTtrialType
for ii=1:NbTtrialType
    ii
    if psth_trial_type{ii} == 'all'
        trigtimes = trigtimes_all ;
    elseif psth_trial_type{ii} == 'cor' % error NOlick
        trigtimes = trigtimes_cor ;
    elseif psth_trial_type{ii} == 'imp'
        trigtimes = trigtimes_imp ;
    elseif psth_trial_type{ii} == 'opt' % error NOlick
        trigtimes = trigtimes_opt ;
    elseif psth_trial_type{ii} == 'cCL'
        trigtimes = trigtimes_cCL ;
    elseif psth_trial_type{ii} == 'cCR'
        trigtimes = trigtimes_cCR ;
    elseif psth_trial_type{ii} == 'iCL'
        trigtimes = trigtimes_iCL ;
    elseif psth_trial_type{ii} == 'iCR'
        trigtimes = trigtimes_iCR ;
    elseif psth_trial_type{ii} == 'oCR' % opto correct Choice right
        trigtimes = trigtimes_oCR ;
    elseif psth_trial_type{ii} == 'oCL'
        trigtimes = trigtimes_oCL ;
    elseif psth_trial_type{ii} == 'ocC' % opto correct Choices l+r
        trigtimes = trigtimes_ocC ;
    elseif psth_trial_type{ii} == 'oNO' % opto NOlick
        trigtimes = trigtimes_oNO ;
    elseif psth_trial_type{ii} == 'ePL' % error Side Puff left
        trigtimes = trigtimes_ePL ;
    elseif psth_trial_type{ii} == 'ePR'
        trigtimes = trigtimes_ePR ;
    elseif psth_trial_type{ii} == 'eNO' % error NOlick
        trigtimes = trigtimes_eNO ;
    end
    
    %% Compute SDF for each trials
    sdf=[];
    nbtrial= max(size(trigtimes));
    for tt=1:nbtrial
        fr=1; %if '1', normalizes to firing rate (Hz); if '0', does nothing (default)
        binsz=1; %ms  bin size of psth (default: 1 ms)
        [psth_tt] = mpsth(spxtimes, trigtimes(tt), 'pre', pre, 'post', post, 'fr',fr, 'binsz', binsz, 'chart',0, 'tb',0);
        
        ftype = 'Gauss' % boxcar, Gauss, exp, exGauss
        %             w = 100 %ms
        w=GaussSmooth;
        [sdf_tt kernel] = msdf(psth_tt,ftype,w);
        sdf=[sdf sdf_tt];
    end
    
    %% PLOT RASTER
    subplot(NbTtrialType+2,1,ii)
    [psth trialspx] = mpsth(spxtimes, trigtimes, 'pre', pre, 'post', post, 'fr',fr, 'binsz', binsz, 'chart',0, 'tb',0);
    rastmat = zeros(numel(trialspx),pre+1+post);
    timevec = -pre:1:post;
    for i = 1:numel(trialspx)
        rastmat(i,trialspx{i}+pre+1) = 1;
        plot(timevec(rastmat(i,:)~=0),rastmat(i,rastmat(i,:)~=0)*i,'.','Color',col{ii},'MarkerSize',5),hold on
        legend({['trial ' psth_trial_type{ii}]}, 'FontSize',11,'FontWeight', 'normal')
    end
    gca, axis([-pre+10 post+10 0.5 numel(trialspx)+0.5]);
    xlabel(['time from ' psth_trig_evt '(ms)']);
    ylabel('trials');
    xlim([-pre+250 post-250])
    ax = gca; ax.Visible = 'off';  % Set the 'visible' property 'off'
    
    %% PLOT SDF
    subplot(NbTtrialType+2,1,[NbTtrialType+1 NbTtrialType+2]),
    sdf_mean= mean(sdf');
    sdf_std = std(sdf');
    sdf_sem = sdf_std/sqrt(nbtrial);
    
    fstr= col{ii};
    hold on, plot([-pre:1:post], sdf_mean,'color',col{ii},'LineWidth',2);
    hold on, plotshaded([-pre:1:post], [sdf_mean+(K*sdf_sem); sdf_mean-(K*sdf_sem)] ,fstr);
    %     hold on, plotshaded([-pre:1:post], [sdf_mean+sdf_std; sdf_mean-sdf_std] ,fstr);
    xlim([-pre+250 post-250])
    hold on, plot([-1500 -1500], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5)
    hold on, plot([-750 -750], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5)
    hold on, plot([0 0], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5)
    hold on, plot([750 750], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5)
    hold on, plot([1500 1500], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5)
    
    ylabel(['fr (Hz)'],'FontSize',11,'FontWeight', 'normal')
    xlabel(['time from ' psth_trig_evt ' start (ms)'], 'FontSize',11,'FontWeight', 'normal');
    
end
title([MouseID ' ' Day ' ' ChanID ' clust#' num2str(CLUST) ' ntrial=' num2str(nbtrial)],'FontSize',11,'FontWeight', 'bold')

tag_trtype = [];
for ptt=1:NbTtrialType
    tag_trtype = [tag_trtype psth_trial_type{ptt}];
end


% %     SAVING
% mkdir('D:\JC_Figures\SUA\SDF\Fig_Article')
% saveas(gcf, ['D:\JC_Figures\SUA\SDF\Fig_Article\' MouseID '_' Day  '_' ChanID '_clust#' num2str(CLUST) '_SDF_' tag_trtype '_' psth_trig_evt ], 'png')
% saveas(gcf, ['D:\JC_Figures\SUA\SDF\Fig_Article\' MouseID '_' Day  '_' ChanID '_clust#' num2str(CLUST) '_SDF_' tag_trtype '_' psth_trig_evt ], 'eps')
%
% mkdir('D:\JC_Figures\Fig_Article')
% saveas(gcf, ['D:\JC_Figures\Fig_Article\' MouseID '_' Day  '_' ChanID '_clust#' num2str(CLUST) '_SDF_' tag_trtype '_' psth_trig_evt ], 'png')
% saveas(gcf, ['D:\JC_Figures\Fig_Article\' MouseID '_' Day  '_' ChanID '_clust#' num2str(CLUST) '_SDF_' tag_trtype '_' psth_trig_evt ], 'eps')


