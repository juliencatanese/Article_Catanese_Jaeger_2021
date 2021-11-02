% plot_RASTER_Spk_JCfun
% Require the use of 'dat2mat_JC_Script.m' to generate .mat files
% Require the use of 'wave_clus.m' to generate .mat files
% By Julien Catanese in JaegerLab 2017
%%
close all;
% dat2mat_JC_Script


%% Define trials of interest (sort by trials types)
% trial_sort = 'all'
trial_sort = 'cle' % clean
% trial_sort = 'lef'
% trial_sort = 'rig'
% trial_sort = 'cor'

%% load data
load('info.mat')
load('evt.mat')
load('time.mat')
sr=info.info_freq_parameters.board_dig_in_sample_rate;


%% define trial idx
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);

idx_trial_start = trig_end(1:end-1) -(1*sr);
idx_trial_end = trig_st(2:end)      -(1*sr);

Ntrials= max(size(idx_trial_start))

%% define events
evt_valve = evt_rwd_L + evt_rwd_R;
evt_lick = evt_lick_L + evt_lick_R;
evt_puff =evt_puff_L + evt_puff_R;
evt_preLick = evt_delay + evt_puff;


%% Trial selection : correct & CLEAN trials (Sporadic Licks before GO cue are removed)
tr_ID_Resp_good = [];  idx_Resp_good= [];  idx_Go_clean = []; idx_Puff_clean = [];

for tr=1:Ntrials   % Select trials with correct resp Lick AND no lick during Puff and delay periods.
    idx_tr_st = idx_trial_start(tr);
    idx_tr_end= idx_trial_end(tr);
    if (sum(evt_valve(idx_tr_st:idx_tr_end) & evt_lick(idx_tr_st:idx_tr_end))>0) & (sum(evt_preLick(idx_tr_st:idx_tr_end) & evt_lick(idx_tr_st:idx_tr_end))== 0)
        tr_ID_Resp_good = [tr_ID_Resp_good tr];
        disp(['Correct Resp at trial ' num2str(tr)]);
        idx_Resp_g = []; idx_go=[]; idx_puff = []; 
        idx_Resp_g = find(diff(evt_valve(idx_tr_st:idx_tr_end) & evt_lick(idx_tr_st:idx_tr_end))>0)+1; %Start resp & valve
        idx_go = find(diff(evt_delay(idx_tr_st:idx_tr_end))<0); % End of delay
        idx_puff = find(diff(evt_puff(idx_tr_st:idx_tr_end))>0)+1; % Start Puff
        
        if max(size(idx_Resp_g))>1
            idx_Resp_g = min(idx_Resp_g);
            
            disp(['Multiple Resp at trial ' num2str(tr)]);
        end
        % idx within TRIAL
        idx_Resp_good = [idx_Resp_good idx_Resp_g];
        idx_Go_clean = [idx_Go_clean idx_go];% idx within TRIAL of the Response Lick that trigger rwd (the initiation)
        idx_Puff_clean = [idx_Puff_clean idx_puff];
    end
    disp(['nothing at trial ' num2str(tr)]);
end
%%
disp(['Number of Clean Correct trials = ' num2str(max(size( tr_ID_Resp_good )))])
disp(['Percent of Clean Correct trials = ' num2str(100*(max(size( tr_ID_Resp_good ))/Ntrials)) '%'])

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
end

%% loop trough Channels and Clusters
% The variable cluster_class has 2 columns and nspk rows ?nspk is the number of spikes?. The
% first column is the cluster class, with integers denoting the clusters membership and a value of 0 for those spikes not assigned
% to any cluster. The second column is the spike times in ms.

clust_file = dir('times_*_sub*.mat');

for ncluf = 1: max(size(clust_file)) %  channel loop
    load(clust_file(ncluf).name); %  e.g.: load('times_S2Ch6_sub.mat')
    for CLUST= 1:max(cluster_class(:,1)) %Cluster loop within each channel
        spktime_ms=cluster_class(find(cluster_class(:,1)==CLUST),2);
        % convert time to sec
        spkt = spktime_ms/10^3;
        % get the idx
        idxspk=round(spkt*sr);
        % Create a spk boolean vec ([0 0 0 1 0 0 1 0 0 1 0 1 ]) same size as time and evt
        spkbool = zeros(size(time));
        spkbool(idxspk)=1;
        
        Nspk=max(size(spkt))
        
        %% transform 0 into NaN
        spkbool(~spkbool)=NaN;
        
        
        
        
        for Nb = 1:Nblock
%             Nb
            if Ntrial < trialperBlock*Nb
                endblock = Ntrial;
            else
                endblock = trialperBlock*Nb;
            end
            
            if Nb<10
                zz = '0';
            else
                zz='';
            end
            
            close all;
            figure,
            pause(0.1)
            for tr=1+(trialperBlock*(Nb-1)):endblock
                pause(0.01)
                time_tr = time(idx_trial_start(tr):idx_trial_end(tr))-time(idx_trial_start(tr));
                
                plot(time_tr,tr*evt_delay(idx_trial_start(tr):idx_trial_end(tr)),'color',[0.5 0.5 0.5],'LineWidth',2)
                xlim([0,4])
                
                hold on,
                plot(time_tr,tr*evt_puff_L(idx_trial_start(tr):idx_trial_end(tr)),'m') %m
                hold on,
                plot(time_tr,tr*evt_puff_R(idx_trial_start(tr):idx_trial_end(tr)),'c') %c
                hold on,
                
                plot(time_tr,tr*evt_rwd_L(idx_trial_start(tr):idx_trial_end(tr)),'*k')
                hold on,
                plot(time_tr,tr*evt_rwd_R(idx_trial_start(tr):idx_trial_end(tr)),'*k')
                hold on,
                
                plot(time_tr,tr*evt_lick_L(idx_trial_start(tr):idx_trial_end(tr)),'.r') %color',[1 0.45 0.65])
                hold on,
                plot(time_tr,tr*evt_lick_R(idx_trial_start(tr):idx_trial_end(tr)),'.b')
                hold on,
                
                plot(time_tr,tr*evt_opto(idx_trial_start(tr):idx_trial_end(tr)),'y','LineWidth',2)
                hold on,
                
                plot(time_tr,tr*spkbool(idx_trial_start(tr):idx_trial_end(tr)),'.k')
                hold on,
                %         plot(time_tr,tr*evt_puff_L(idx_trial_start(tr):idx_trial_end(tr)),'r','LineWidth',2) %m
                %         hold on,
                
                xlabel('time (s)')
                ylabel('# trials')
                
                title( [info.info_notes.MouseID ' '  info.info_notes.Day(1:end-3) ' ' clust_file(ncluf).name(7:12) ' CLUST#' num2str(CLUST) ' block#' zz num2str(Nb) ])
                %         xlim([0 4.2])
                %         plot(Spktime, ones(size(Spktime)),'*r')
                
                
            end
            %     xlim([0,13])
            saveas(gcf,[ info.info_notes.MouseID '_'  info.info_notes.Day(1:end-4)  '_SpkRASTER_' clust_file(ncluf).name(7:11) '_CLUST#' num2str(CLUST) '_block#' zz num2str(Nb) '_' trial_sort 'trial_v2' ] ,'tif')
            close gcf
        end
    end
end

disp('COMPLETE SUCCESSFULLY: SEE FIGURES IN FOLDER')
% %% define Licks idx
% idx_lick_L = find(diff(evt_lick_L(idx_trial_start:idx_trial_end))>0);
% idx_Lick_R = find(diff(evt_lick_R)>0);

