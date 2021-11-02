% Pop_Seq_DistMatrix_1session_GOcent_corr_mlibJCscript.m
% plot normalized SDF for each neurons in the session in each row 
% color code is the normalized firing rate 
% Xaxis is time in ms centered on Gocue  
% written by Julien Catanese   % use mlib library 
% last updated JC 11/21/2018 


%% Define PSTH parameter
try
    psth_trig_evt =  Center
    psth_trial_type = {psth_trial_type};
    NbTtrialType = max(size(psth_trial_type))
    col={color}
catch
    psth_trig_evt =  'Delay'; % center for SDF
    psth_trial_type = {'cor'};
    NbTtrialType = max(size(psth_trial_type));
    col={'k'}
    pre=1500;
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
for tr=1:Ntrials
    time_tr = time(idx_trial_start(tr):idx_trial_end(tr));
    if psth_trig_evt=='Delay'
        delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
        idx_delay_st = find(diff(delay_tr)>0)  % start of the delay
        time_delay_st = time_tr(idx_delay_st);
        trigtimes = [trigtimes time_delay_st]; % in sec
    elseif psth_trig_evt=='GoCue'
        delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
        idx_GO_st = find(diff(delay_tr)<0); % end of the delay
        time_GO_st = time_tr(idx_GO_st);
        trigtimes = [trigtimes time_GO_st]; % in sec
    else
        disp(' WRONG INPUT ARG pick one: psth_trig_evt = GoCue, APuff, Delay, Licks or Valve' )
    end
end

%%
trigtimes_all=trigtimes';
trigtimes_cCL = trigtimes(trial.idx_correct_L);
trigtimes_cCR = trigtimes(trial.idx_correct_R);
trigtimes_cor = sort([trigtimes_cCR trigtimes_cCL]);

trigtimes = trigtimes_cor ;
%% GET spxtimes = SPIKE times in Sec (1 vector of times)
close all,
clust_file = dir('times_*S*Ch*_sub.mat');

psthall = [];
sdfall = []; 
for ncluf = 1:max(size(clust_file)) %  channel loop
    
    load(clust_file(ncluf).name); %  e.g.: load('times_S2Ch6_man.mat')
    Nclust = max(cluster_class(:,1));
    
 
    for CLUST= 1:Nclust %Cluster loop within each channel
        idx_spk = find(cluster_class(:,1)==CLUST);
        spxtimes_ms = cluster_class(idx_spk,2);
        spxtimes = spxtimes_ms/10^3; % convert time to sec
% 
%         pre=1500; %ms
%         post=1500;%ms
        fr=1; %if '1', normalizes to firing rate (Hz); if '0', does nothing (default)
        binsz=1; %ms  bin size of psth (default: 1 ms)
        [psth trialspx] = mpsth(spxtimes, trigtimes, 'pre', pre, 'post', post, 'fr',fr, 'binsz', binsz, 'chart',0);
        
        figure (1),
        ftype = 'Gauss' % boxcar, Gauss, exp, exGauss
        w = 100 %ms
        [sdf kernel] = msdf(psth,ftype,w);
        hold on, plot(sdf(:,1), sdf(:,2),'color','b','LineWidth',2 );
        
       
        psthall = [psthall; psth(:,2)'];
        sdfall = [sdfall; sdf(:,2)'];
    end
end

%%
% DD=sdfall(1,:)./max(sdfall(1,:))

AA = max(sdfall');
NN = (sdfall'./AA);
% figure, plot(NN)
for irow= 1:size(NN,2)
MM(irow) =  min(find(NN(:,irow)==1))
end
BB = [MM; NN];
sort_sdfall = sortrows(BB','descend');
figure, imagesc(sort_sdfall(:,2:end))
colorbar, 
caxis([0.5 1])
title([MouseID ' ' Day ])



%% saving
saveas(gcf, [MouseID '_' Day '_Pop_SequenceMatriceSorted.png'])
saveas(gcf, ['D:\JC_Figures\pop\' MouseID '_' Day '_Pop_SequenceMatriceSorted.png'])


% ALL_SDFALL = [ALL_SDFALL; sort_sdfall] 
clear AA NN MM BB sortsdfall



