% Stats_plot_LvR_epoch_JC_Script
% Plot Histogram for STATS Firing rate Left vs Right trials
% This function use the mlib toolbox for mraster, mpsth and msdf.
% Written by Julien Catanese, 10/03/2018, in JaegerLab.
% Last Updated: 10/08/2018
psth_center_evtID = 'Delay'
% Other_trial_type = []
clearvars -except  psth_center_evtID Other_trial_type FolderID MouseID
close all,

%% Define PSTH parameter

% if nargin==0;
%         display('NOT ENOUGH INPUT ARGUMENTS: min=1')
% elseif nargin==1;
disp(['psth Center = ' psth_center_evtID]);
disp('NO Other trial type, plotting LvR only');
trial_type_list = {'cCL', 'cCR'};
Special_trial= 'x';
% elseif nargin==2;
%      display(['Other trial type = ' Other_trial_type]);
%      trial_type_list = {'cCL', 'cCR', Other_trial_type};
%      Special_trial= Other_trial_type;
% elseif nargin>2;
%     display('TOO MANY INPUT ARGUMENTS: max=1');
% end

psth_trig_evt =  psth_center_evtID; % 'Delay', % 'APuff' % 'GoCue' % 'Valve' %'Licks'
NbTtrialType = max(size(trial_type_list));
%% load data
load('info.mat');
load('evt.mat');
load('time.mat');
load('Ntrial_type.mat');
sr=info.info_freq_parameters.board_dig_in_sample_rate;

%% GET trigtimes : times of the psth_trig_evt is sec   (1 vector of times)
% define trial idx
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);

idx_trial_start = trig_end(1:end-1); %-(1*sr);
idx_trial_end = trig_st(2:end) ;     %-(1*sr);

Ntrials=trial.Ntrial; 

% define common events
evt_valve = evt_rwd_L + evt_rwd_R;
evt_lick = evt_lick_L + evt_lick_R;
evt_puff =evt_puff_L + evt_puff_R;
evt_prelick = evt_delay + evt_puff;


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
        idx_delay_st = find(diff(delay_tr)>0) ; % start of the delay
        time_delay_st = time_tr(idx_delay_st);
        trigtimes = [trigtimes time_delay_st]; % in sec
        
    elseif psth_trig_evt=='GoCue'
        delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
        idx_GO_st = find(diff(delay_tr)<0) ; % end of the delay
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
        
    elseif psth_trig_evt=='Valve'
        valve_tr= evt_valve(idx_trial_start(tr):idx_trial_end(tr));
        idx_valve_st = find(diff(valve_tr)>0) ; % start of the valve
        time_valve_st = time_tr(idx_valve_st);
        trigtimes = [trigtimes time_valve_st]; % in sec
    else
        disp(' WRONG INPUT ARG pick one: psth_trig_evt = APuff, Delay, Licks, Gocue or Valve' )
        
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

%% GET spxtimes = SPIKE times in Sec (1 vector of times)
close all,
clust_file = dir('times_*S*Ch*_sub.mat');
for ncluf = 1:max(size(clust_file)) %  channel loop
    
    load(clust_file(ncluf).name); %  e.g.: load('times_S2Ch6_man.mat')
    Nclust = max(cluster_class(:,1));
    
    for CLUST= 1:Nclust %Cluster loop within each channel
        idx_spk = find(cluster_class(:,1)==CLUST);
        spxtimes_ms = cluster_class(idx_spk,2);
        % convert time to sec
        spxtimes = spxtimes_ms/10^3;
        Nspk=max(size(spxtimes));
        
        nspx_delay_tr = [];
        nspx_sample_tr = [];
        nspx_resp_tr = [];
        
        figure, hold on,
        for ii=1:NbTtrialType
            
            if trial_type_list{ii} == 'all'
                trigtimes = trigtimes_all ;
            elseif trial_type_list{ii} == 'cCL'
                trigtimes = trigtimes_cCL ; %disp('Correct Left trials');
            elseif trial_type_list{ii} == 'cCR'
                trigtimes = trigtimes_cCR ; %disp('Correct Right trials');
            elseif trial_type_list{ii} == 'iCL'
                trigtimes = trigtimes_iCL ;
            elseif trial_type_list{ii} == 'iCR'
                trigtimes = trigtimes_iCR ;
            elseif trial_type_list{ii} == 'oCR' % opto correct Choice right
                trigtimes = trigtimes_oCR ;
            elseif trial_type_list{ii} == 'oCL'
                trigtimes = trigtimes_oCL ;
            elseif trial_type_list{ii} == 'ocC' % opto correct Choices l+r
                trigtimes = trigtimes_ocC ;
            elseif trial_type_list{ii} == 'oNO' % opto NOlick
                trigtimes = trigtimes_oNO ;
            elseif trial_type_list{ii} == 'ePL' % error Side Puff left
                trigtimes = trigtimes_ePL ;
            elseif trial_type_list{ii} == 'ePR'
                trigtimes = trigtimes_ePR ;
            elseif trial_type_list{ii} == 'eNO' % error NOlick
                trigtimes = trigtimes_eNO ;
            end
            
            pre=2000; %ms
            post=2000;%ms
            fr=1; %if '1', normalizes to firing rate (Hz); if '0', does nothing (default)
            binsz=1; %ms  bin size of psth (default: 1 ms)
            
            
            %SAMPLE:  Count Spikes within Puff Sample Epoch for each trial ([delaystart-750] to [delaystart+0] )
            nspx = mnspx(spxtimes,trigtimes,750,0);
            subplot(3,2,ii),
            if sum(nspx)~=0
                hist(nspx,0.5:max(nspx)+0.5); xlim([0 10 ])
                title (['Sample Epoch '  trial_type_list{ii} ' trials' ])
                xlabel('#spk/trial')
                ylabel('#trials')
                if ii==1
                    nspx_sample_left_tr = nspx;
                elseif ii==2
                    nspx_sample_right_tr = nspx;
                end
            end
            
            % DELAY : Count Spikes within Delay Epoch for each trial ([delaystart+0] to [delaystart+750ms] )
            nspx = mnspx(spxtimes,trigtimes,0,750);
            subplot(3,2,ii+2),
            if sum(nspx)~=0
                hist(nspx,0.5:max(nspx)+0.5); xlim([0 10 ]);
                title (['Delay Epoch '  trial_type_list{ii} ' trials' ])
                xlabel('#spk/trial')
                ylabel('#trials')
                if ii==1
                    nspx_delay_left_tr = nspx;
                elseif ii==2
                    nspx_delay_right_tr = nspx;
                end
            end
            
            %RESPONSE:  Count Spikes within Response Epoch for each trial ([delaystart+750] to [delaystart+2000] )
            nspx = mnspx(spxtimes,trigtimes,-750,2000);
            subplot(3,2,ii+4),
            if sum(nspx)~=0
                histogram(nspx,0.5:max(nspx)+0.5, 'Normalization', 'probability')    ;
                title (['Response Epoch '  trial_type_list{ii} ' trials'  ])
                xlabel('#spk/trial')
                ylabel('#trials')
                if ii==1
                    nspx_resp_left_tr = nspx;
                elseif ii==2
                    nspx_resp_right_tr = nspx;
                end
            end
        end
        
        %% STAT: T-TEST2
        Mouse=info.info_notes.MouseID;
        Day= info.info_notes.Day;
        % SAMPLE
        [H,P] = ttest2(nspx_sample_left_tr, nspx_sample_right_tr);
        if H==1
            disp([ Mouse ' ' Day ' ' clust_file(ncluf).name(7:11) ' CLUST#' num2str(CLUST)]);
            disp('SAMPLE Epoch: Significant L-R firing rate ');
            P
            figure,
            %             yyaxis left
            histogram(nspx_sample_left_tr,0.5:max(nspx_sample_left_tr)+0.5, 'Normalization', 'count')
            %             ylabel('Proba')
            hold on,
            %             yyaxis right
            histogram(nspx_sample_right_tr,0.5:max(nspx_sample_right_tr)+0.5, 'Normalization', 'count')
            title ([Mouse ' ' Day '   ' clust_file(ncluf).name(7:11) ' CLUST#' num2str(CLUST) '   SAMPLE  P=' num2str(P)])
            ylabel('Count (#trials)')
            xlabel('#spk/trial')
            legend ('Ipsi sample', 'Contra sample')
            saveas(gcf, ['fig2print_' clust_file(ncluf).name(7:11) ' clust#' num2str(CLUST) '_ttest_LvR_sample'],'png')
        end
        
        % DELAY
        [H,P] = ttest2(nspx_delay_left_tr, nspx_delay_right_tr);
        if H==1
            disp([ Mouse ' ' Day ' ' clust_file(ncluf).name(7:11) ' CLUST#' num2str(CLUST)]);
            disp('DELAY Epoch: Significant L-R firing rate ');
            P
            figure, hold on,
            histogram(nspx_delay_left_tr,0.5:max(nspx_delay_left_tr)+0.5, 'Normalization', 'count')
            histogram(nspx_delay_right_tr,0.5:max(nspx_delay_right_tr)+0.5, 'Normalization', 'count')
            title ([Mouse ' ' Day '   ' clust_file(ncluf).name(7:11) ' CLUST#' num2str(CLUST) '   DELAY  P=' num2str(P)])
            ylabel('Count (#trials)')
            xlabel('#spk/trial')
            legend ('Ipsi delay', 'Contra delay')
            saveas(gcf, ['fig2print_' clust_file(ncluf).name(7:11) ' clust#' num2str(CLUST) '_ttest_LvR_delay'],'png')
        end
        
        % RESP
        [H,P] = ttest2(nspx_resp_left_tr, nspx_resp_right_tr);
        if H==1
            disp([ Mouse ' ' Day ' ' clust_file(ncluf).name(7:11) ' CLUST#' num2str(CLUST)]);
            disp('RESP Epoch: Significant L-R firing rate ');
            P
            figure, hold on,
            histogram(nspx_resp_left_tr,0.5:max(nspx_resp_left_tr)+0.5, 'Normalization', 'count')
            histogram(nspx_resp_right_tr,0.5:max(nspx_resp_right_tr)+0.5, 'Normalization', 'count')
            title ([Mouse ' ' Day '   ' clust_file(ncluf).name(7:11) ' CLUST#' num2str(CLUST) '   RESPONSE  P=' num2str(P)])
            ylabel('Count (#trials)')
            xlabel('#spk/trial')
            legend ('Ipsi response', 'Contra response')
            saveas(gcf, ['fig2print_' clust_file(ncluf).name(7:11) ' clust#' num2str(CLUST) '_ttest_LvR_response'],'png')
        end
        
        
    end
end

