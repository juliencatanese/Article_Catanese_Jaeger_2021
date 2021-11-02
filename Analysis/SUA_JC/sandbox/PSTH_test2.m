% PSTH_test2
% Require the use of 'JC_AnalysisDJlab_Master_Script.m' the SXChX_sub.mat files 
% Require the use of 'Behavior_Master_Script.m' the Ntrials_type.mat files
% Require the use of 'wave_clus.m' to generates the times_SXChX_sub.mat files  
% Require to complete the manual Spike sorting. 
% By Julien Catanese 10/03/2018

%%
clear all, close all;

%% Define trials of interest (sort by trials types)

% trial_sort = 'cle'
% trial_sort = 'all'
trial_sort = 'lef'
% trial_sort = 'rig'
% trial_sort = 'cor'

center = 'puff'
% center = 'GOcu'
% center = 'resp'

Nsec = 3 % in second around the center
% Nsec = 2
% Nsec=1



%% load data
load('info.mat');
load('evt.mat');
load('time.mat');
load('Ntrial_type.mat');

sr=info.info_freq_parameters.board_dig_in_sample_rate;

%% define trial idx
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);

idx_trial_start = trig_end(1:end-1) -(1*sr);
idx_trial_end = trig_st(2:end)      -(1*sr);

Ntrials= max(size(idx_trial_start)); if Ntrials==trial.Ntrial; disp('good Ntrial'); else forcestop; end ;

%% define events
evt_valve = evt_rwd_L + evt_rwd_R;
evt_lick = evt_lick_L + evt_lick_R;
evt_puff =evt_puff_L + evt_puff_R;
evt_preLick = evt_delay + evt_puff;


%% Trial selection : correct & CLEAN trials (Sporadic Licks before GO cue are removed)
tr_ID_Resp_good = [];  idx_Resp_clean= [];  idx_Go_clean = []; idx_Puff_clean = [];

for tr=1:75-1   % Select trials with correct resp Lick AND no lick during Puff and delay periods.
    idx_tr_st = idx_trial_start(tr);
    idx_tr_end= idx_trial_end(tr);
    %     select only correct trial with no licks before Go cue
    if (sum(evt_valve(idx_tr_st:idx_tr_end) & evt_lick(idx_tr_st:idx_tr_end))>0) & (sum(evt_preLick(idx_tr_st:idx_tr_end) & evt_lick(idx_tr_st:idx_tr_end))== 0);
        
        disp(['Correct Resp at trial ' num2str(tr)]);
        
        tr_ID_Resp_good = [tr_ID_Resp_good tr];
        
        % idx within TRIAL
        idx_Resp = min(find(diff(evt_valve(idx_tr_st:idx_tr_end) & evt_lick(idx_tr_st:idx_tr_end))>0)+1); %Start resp & valve
        idx_go = find(diff(evt_delay(idx_tr_st:idx_tr_end))<0); % End of delay
        idx_puff = find(diff(evt_puff(idx_tr_st:idx_tr_end))>0)+1; % Start Puff
        
        % idx within TRIAL
        idx_Resp_clean = [idx_Resp_clean idx_Resp];
        idx_Go_clean = [idx_Go_clean idx_go];% idx within TRIAL of the Response Lick that trigger rwd (the initiation)
        idx_Puff_clean = [idx_Puff_clean idx_puff];
    end
    disp(['nothing at trial ' num2str(tr)]);
end
%%
disp(['Number of Clean Correct trials = ' num2str(max(size( tr_ID_Resp_good )))]);
disp(['Percent of Clean Correct trials = ' num2str(100*(max(size( tr_ID_Resp_good ))/Ntrials)) '%']);

%% prepare plot
% figure, hold on,
% Replace zeros by NaN (for further figures)
evt_delay(~evt_delay)=NaN;
evt_opto(~evt_opto)=NaN;

evt_lick_L(~evt_lick_L)=NaN;
evt_lick_R(~evt_lick_R)=NaN;

evt_puff_L(~evt_puff_L)=NaN;
evt_puff_R(~evt_puff_R)=NaN;

evt_rwd_L(~evt_rwd_L)=NaN;
evt_rwd_R(~evt_rwd_R)=NaN;

%% define #trials per Block (per figure)

trialperBlock = 70

if trial_sort == 'all'
    Ntrial = max(size(idx_trial_start));
    Nblock = ceil(Ntrial/trialperBlock);
elseif trial_sort == 'cle'
    Ntrial= max(size(tr_ID_Resp_good));
    Nblock = ceil(Ntrial/trialperBlock);
    idx_trial_start=idx_trial_start(tr_ID_Resp_good);
    idx_trial_end=idx_trial_end(tr_ID_Resp_good);
elseif trial_sort == 'lef'
    Ntrial= max(size(trial.idx_correct_L));
    Nblock = ceil(Ntrial/trialperBlock);
    idx_trial_start=idx_trial_start(trial.idx_correct_L);
    idx_trial_end=idx_trial_end(trial.idx_correct_L);
elseif trial_sort == 'rig'
    Ntrial= max(size(trial.idx_correct_R));
    Nblock = ceil(Ntrial/trialperBlock);
    idx_trial_start=idx_trial_start(trial.idx_correct_R);
    idx_trial_end=idx_trial_end(trial.idx_correct_R);
end


%% loop trough Channels and Clusters
% The variable cluster_class has 2 columns and nspk rows ?nspk is the number of spikes?. The
% first column is the cluster class, with integers denoting the clusters membership and a value of 0 for those spikes not assigned
% to any cluster. The second column is the spike times in ms.

clust_file = dir('times_*Ch*_sub.mat');
for ncluf = 1: max(size(clust_file)) %  channel loop
    load(clust_file(ncluf).name); %  e.g.: load('times_S2Ch6_man.mat')
    for CLUST= 1:max(cluster_class(:,1)) %Cluster loop within each channel
        spktime_ms=cluster_class(find(cluster_class(:,1)==CLUST),2);
        % convert time to sec
        spkt = spktime_ms/10^3;
        % get the idx
        idxspk=round(spkt*sr);
        % Create a spk boolean vec ([0 0 0 1 0 0 1 0 0 1 0 1 ]) same size as time and evt
        spkbool = zeros(size(time));
        spkbool(idxspk)=1;
        
        Nspk=max(size(spkt));
        
        %% transform 0 into NaN
        %         spkbool(~spkbool)=NaN;
        
        for tr=1:Ntrial
            if center == 'puff';
                idx_center = idx_trial_start(tr) + idx_Puff_clean(tr);
            elseif center == 'resp';
                idx_center = idx_trial_start(tr) + idx_Resp_clean(tr);
            elseif center == 'Gocu';
                idx_center = idx_trial_start(tr) + idx_GO_clean(tr);
            end
            
            data_center(tr,:)=spkbool(idx_center - (Nsec*sr) : idx_center + (Nsec*sr));
            time_center(tr,:)=time(idx_center - (Nsec*sr) : idx_center + (Nsec*sr));
            evt_puff_center(tr,:)=evt_puff(idx_center - (Nsec*sr) : idx_center + (Nsec*sr));
            evt_delay_center(tr,:)=evt_delay(idx_center - (Nsec*sr) : idx_center + (Nsec*sr));
            evt_lick_center(tr,:)=evt_lick(idx_center - (Nsec*sr) : idx_center + (Nsec*sr));
        end
    
      
    %% Defines Bins for PSTH %%%
    %%
    timebin = 0.25  % second
    nbin = size(data_center,2)/sr/timebin;
    binsize = size(data_center,2)/nbin;
    bin_mat =  ones(size(data_center));
  %%  create the bin_mat [1 1 1 2 2 2 3 3 3 ... 24 24 24 ; 1 1 1 2 2 2 3 3 3 ... 24 24 24 ; ...] x 30 trials (24 bins per row, 30 trials, 1 trial per row) 
    for nb = 1:nbin
        bin_mat(:,(1+((nb-1)*binsize)):(binsize+((nb-1)*binsize))) = nb;  
    end
   disp('bin_mat done');

%% cout the number of spikes in each bin 
      DMb = data_center .* bin_mat;
      ePc = evt_puff_center(1,:).*bin_mat(1,:);
      eDc=  evt_delay_center(1,:).*bin_mat(1,:);
      eLc=  evt_lick_center(1,:).*bin_mat(1,:);
      
      binpuff=[], bindelay=[], binlick=[], bincount=[];  
    for nb = 1:nbin
           bincount(nb) =  size(find(DMb==nb),1); 
           binpuff(nb) =  size(find(ePc==nb),2)>1; 
           bindelay(nb)= size(find(eDc==nb),2)>1; 
           binlick(nb)= size(find(eLc==nb),2)>1;
    end
    
    binpuff(~binpuff) = NaN  ;
    bindelay(~bindelay) = NaN  ;
    binlick(~binlick) = NaN;
    
    figure, hold on,
    bar(bincount);
    
    plot(binpuff*1,'r');
    plot(bindelay*1,'g');
    plot(binlick*1,'b');
    legend ('spk', 'puff','delay','licks');
    xlabel ('bins (0.25s)');
    ylabel ('#spikes')  ;
    title(['PSTH centered on Puff onset ' clust_file(ncluf).name(7:12) ' CLUST#' num2str(CLUST)]);
    
    saveas(gcf,[ info.info_notes.MouseID '_'  info.info_notes.Day(1:end-4)  '_PSTH_' clust_file(ncluf).name(7:11) '_CLUST#' num2str(CLUST) '' '_center' center '_' trial_sort 'trial' ] ,'tif');

    end
    
end
 

disp('COMPLETE SUCCESSFULLY: SEE FIGURES IN FOLDER')
