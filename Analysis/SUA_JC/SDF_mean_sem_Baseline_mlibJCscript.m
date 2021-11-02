% SDF_mean_sem_Baseline_mlibJCscript
% plot normalized mean SDF with SEM obtain by computing each SDF per trial 
% Default parameters: 
%     psth_trig_evt =  'Delay'; % psth centered on delay start 
%     psth_trial_type = {'cor'}; % only correct trials
%     pre=2250+750; % 2250ms baseline + 750ms air puff
%     post=1500;  % 750ms delay + 750ms response
% dependency: mlib6 library
% written by Julien Catanese 11/27/2018 
% last updated JC 11/28/2018


%% Define PSTH parameter
try
    psth_trig_evt
    pre
    post
    psth_trial_type = {psth_trial_type};
    NbTtrialType = max(size(psth_trial_type));
catch % Defaults
    disp('Using default parameters: center Delay; correct trials; pre=post=1500ms');
    psth_trig_evt =  'Delay'; % center for SDF
    psth_trial_type = {'cor'};
    NbTtrialType = max(size(psth_trial_type));
    pre=2250+750;
    post=1500;
end


%% load data
load('info.mat');  MouseID= info.info_notes.MouseID; Day = info.info_notes.Day;
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
Ntrials= trial.Ntrial;

%% define trigtimes for each trials
trigtimes=[];
for tr=1:Ntrials;
    time_tr = time(idx_trial_start(tr):idx_trial_end(tr));
    if psth_trig_evt=='Delay';
        delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
        idx_delay_st = find(diff(delay_tr)>0);  % start of the delay
        time_delay_st = time_tr(idx_delay_st);
        trigtimes = [trigtimes time_delay_st]; % in sec
    elseif psth_trig_evt=='GoCue';
        delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
        idx_GO_st = find(diff(delay_tr)<0); % end of the delay
        time_GO_st = time_tr(idx_GO_st);
        trigtimes = [trigtimes time_GO_st]; % in sec
    else
        disp(' WRONG INPUT ARG pick one: psth_trig_evt = GoCue, APuff, Delay, Licks or Valve' );
    end
end

%% Get final trigtimes
trigtimes_all=trigtimes';
trigtimes_cCL = trigtimes(trial.idx_correct_L);
trigtimes_cCR = trigtimes(trial.idx_correct_R);
trigtimes_cor = sort([trigtimes_cCR trigtimes_cCL]);

trigtimes = trigtimes_cor ;

%% GET spxtimes = SPIKE times in Sec (1 vector of times)
clust_file = dir('times_*S*Ch*_sub.mat');
sdf_norm_mean_ALL = []; sdf_norm_sem_ALL = [];
nchan=max(size(clust_file));
for ncluf = 1:nchan %  channel loop
    
    load(clust_file(ncluf).name); %  e.g.: load('times_S2Ch6_man.mat')
    Nclust = max(cluster_class(:,1));
    
    %     close all;
    for CLUST= 1:Nclust %Cluster loop within each channel
        idx_spk = find(cluster_class(:,1)==CLUST);
        spxtimes_ms = cluster_class(idx_spk,2);
        spxtimes = spxtimes_ms/10^3; % convert time to sec
        if figplot==1; figure(CLUST+1), close(CLUST+1); end; 
        
        %% SDF for each trial
        psth=[]; sdf=[];
        nbtrial = max(size(trigtimes));
        for tt=1:nbtrial
            
            fr=1; %if '1', normalizes to firing rate (Hz); if '0', does nothing (default)
            binsz=1; %ms  bin size of psth (default: 1 ms)
            [psth_tt trialspx] = mpsth(spxtimes, trigtimes(tt), 'pre', pre, 'post', post, 'fr',fr, 'binsz', binsz, 'chart',0, 'tb',0);
            
            ftype = 'Gauss' % boxcar, Gauss, exp, exGauss
            w = 100 %ms
            [sdf_tt kernel] = msdf(psth_tt,ftype,w);
            if figplot==1;
                figure(CLUST+1);
                hold on, plot([-pre:1:post],sdf_tt,'color','b','LineWidth',1 );
            end
            sdf=[sdf sdf_tt];
        end
        %%
        sdf=sdf';
        sdf_mean= mean(sdf);
        sdf_std= std(sdf);
        sdf_sem = std(sdf)/sqrt(nbtrial);
        if figplot==1
            hold on, plot([-pre:1:post],sdf_mean,'color','r','LineWidth',2 );
            hold on, plot([-pre:1:post],sdf_mean+sdf_sem,'--m','LineWidth',2 );
            hold on, plot([-pre:1:post],sdf_mean-sdf_sem,'--m','LineWidth',2 );
        end
        %% normalize to the pic
        sdf_norm_mean= sdf_mean./max(sdf_mean);
        sdf_norm_std= sdf_std./max(sdf_std);
        sdf_norm_sem = sdf_norm_std/sqrt(nbtrial);
        
        if figplot==1
            figure,
            hold on, plot([-pre:1:post],sdf_norm_mean,'r','LineWidth',2);
            hold on, plot([-pre:1:post],sdf_norm_mean+sdf_norm_sem,'--m','LineWidth',2 );
            hold on, plot([-pre:1:post],sdf_norm_mean-sdf_norm_sem,'--m','LineWidth',2 );
        end
        
        sdf_norm_mean_ALL = [sdf_norm_mean_ALL ; sdf_norm_mean];
        sdf_norm_sem_ALL = [sdf_norm_sem_ALL ; sdf_norm_sem];
    end
end
%%
clear SNMA SNSA IDX sort_sdfall

SNMA = sdf_norm_mean_ALL;
SNSA = sdf_norm_sem_ALL;
for irow= 1:size(SNMA,1)
    IDX(irow) =  min(find(SNMA(irow,:)==1));
end
SNMA = [IDX' SNMA];
SNSA = [IDX' SNSA];

sort_sdfall = sortrows(SNMA,'ascend');
sort_sdfall_sem = sortrows(SNSA,'ascend');

figure, imagesc(sort_sdfall(:,2:end));
colorbar,;
caxis([0.5 1]);
title([MouseID ' ' Day ]);



