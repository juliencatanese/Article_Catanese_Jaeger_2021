% Pop_Seq_DistMatrix_1session_GOcent_corr_mlibJCscript.m
% Plot the pop sequence matrice for 1 session
% uses SDF normalized to their peack for each neurons (each row of the Matrice)
% Xaxis is time in ms  center on psth_trig_evt
% color code is the normalized firing rate
%  DEFAULT PARAM:
%     psth_trig_evt =  'Delay'; % center for SDF
%     psth_trial_type = [{'cCR'} {'cCL'}]
%     pre=1500;
%     post=1500;
% dependency : mlib6 library
% written by: Julien Catanese
% last updated JC 11/28/2018



% function SDF_Raster_1cell_mlibJCfun(Sh, Ch, CLUST, psth_trig_evt , psth_trial_type, col, pre, post)
% Default INPUT Paramnters
%     psth_trig_evt =  'Delay'; % psth center on Delay start; can also use 'GoCue'
%     psth_trial_type = {'cor'}; % All correct trials; to spearate Left Right: {'cCL', 'cCR'} for Left Right
%     col={'k'} ; or col={'r', 'b'}
%     pre=1500;
%     post=1500;
%  by JC 11/22/2018

%% Define NbTtrialType
try
    pre
    post
    col
    psth_trig_evt
    psth_trial_type
    NbTtrialType = max(size(psth_trial_type))
    disp('Success input')
catch
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

ncell=0;

%% GET spxtimes = SPIKE times in Sec (1 vector of times)
clust_file = dir('times_*S*Ch*_sub.mat');
sdf_mean_ALL1=[]; sdf_mean_ALL2=[]; sdf_sem_ALL1=[]; sdf_sem_ALL2=[];
ncluf=max(size(clust_file));
for nchan = 1:ncluf %  channel loop
    load(clust_file(nchan).name); %  e.g.: load('times_S2Ch6_man.mat')
    Nclust = max(cluster_class(:,1));
    ncell=ncell+Nclust;
    
    for CLUST= 1:Nclust %Cluster loop within each channel
        idx_spk = find(cluster_class(:,1)==CLUST);
        spxtimes_ms = cluster_class(idx_spk,2);
        spxtimes = spxtimes_ms/10^3; % convert time to sec
        
        for ii=1:NbTtrialType; ii
            
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
                w = 100 %ms
                [sdf_tt kernel] = msdf(psth_tt,ftype,w);
                sdf=[sdf sdf_tt];
            end
            
            if ii==1
                sdf1=sdf'
                sdf_mean1= mean(sdf');
                sdf_sem1 = std(sdf')/sqrt(nbtrial);psth_trig_evt
            elseif ii==2
                sdf2=sdf'
                sdf_mean2= mean(sdf');
                sdf_sem2 = std(sdf')/sqrt(nbtrial);
            end
            
        end
        
        if NbTtrialType == 1
            sdf_mean_ALL1 = [sdf_mean_ALL1 ; sdf_mean1];
            sdf_sem_ALL1 = [sdf_sem_ALL1 ; sdf_sem1];
        elseif NbTtrialType == 2
            sdf_mean_ALL1 = [sdf_mean_ALL1 ; sdf_mean1];
            sdf_sem_ALL1 = [sdf_sem_ALL1 ; sdf_sem1];
            sdf_mean_ALL2 = [sdf_mean_ALL2 ; sdf_mean2];
            sdf_sem_ALL2 = [sdf_sem_ALL2 ; sdf_sem2];
        end
        
    end
end

%%
SMA = []; SSA = []; SNMA =[];  SNSA = []; IDX=[]; sort_sdfall=[]; 

% figure, plot([-pre:1:post], sdf_mean_ALL1-sdf_mean_ALL2);
% 1=contra  2=ipsi
SMA = sdf_mean_ALL1-sdf_mean_ALL2; 
% figure, plot([-pre:1:post], SMA);
SMA=SMA';
Peak = max(abs(SMA)); % negative or positive
SNMA = SMA./Peak;
figure, plot([-pre:1:post], SNMA);
SNMA = SNMA';

for irow= 1:size(SNMA,1)
    IDX(irow) =  min(find(abs(SNMA(irow,:))==1));
    IDX(irow) =  min(find(abs(SNMA(irow,:))==1));
end
SNMA = [IDX' SNMA];
% SNSA = [IDX' SNSA];

sort_sdfall = sortrows(SNMA,'ascend');
% sort_sdfall_sem = sortrows(SNSA,'ascend');

figure, imagesc(sort_sdfall(:,2:end));
colormap('jet')
colorbar,
caxis([-1 1]);
title([MouseID ' ' Day ' contra-ipsi']);

    
    
    
