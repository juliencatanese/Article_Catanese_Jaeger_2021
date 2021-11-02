% plot_Spectro_Mean_Norm_Aligned_JC_Script
% .............
% Spectro can be Aligned on RespLick or Gocue or PuffStart  
% Give 3 aligned visualisations: 1-Normalized Spectro, 2-RawLFPTraces, and 3-Raw Spectro 
% Give the Averaged over correct and clean trials% 
% Give also a trial by trial analysis
% .........
% written by JC May 2018 in Jaeger lab 
% last updated JC 5-30-18
%...........

clear all;
close all;
%% load data
load('info.mat')
load('evt.mat')
load('time.mat')
load('S1Ch6_sub.mat') % data = microV and sr = sampling rate

%% define trial idx
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);

idx_trial_start = trig_end(1:end-1) -(1*sr);
idx_trial_end = trig_st(2:end)      -(1*sr);

Ntrials= max(size(idx_trial_start))

figure, hold on,

%% define the GO epoch and create data_epoch_GO
evt_valve = evt_rwd_L + evt_rwd_R;
evt_lick = evt_lick_L + evt_lick_R;
evt_puff =evt_puff_L + evt_puff_R;
evt_preLick = evt_delay + evt_puff;

clear evt_lick_L
clear evt_lick_R
clear evt_puff_L
clear evt_puff_R
clear evt_rwd_L
clear evt_rwd_R
clear evt_opto

%% Trial selection : correct & CLEAN trials (Sporadic Licks before GO cue are removed)
tr_ID_Resp_good = [];  idx_Resp_good= [];  idx_Go_clean = []; idx_Puff_clean = [];

for tr=1:75-1   % Select trials with correct resp Lick AND no lick during Puff and delay periods.
    idx_tr_st = idx_trial_start(tr);
    idx_tr_end= idx_trial_end(tr);
    if (sum(evt_valve(idx_tr_st:idx_tr_end) & evt_lick(idx_tr_st:idx_tr_end))>0) & (sum(evt_preLick(idx_tr_st:idx_tr_end) & evt_lick(idx_tr_st:idx_tr_end))== 0)
        tr_ID_Resp_good = [tr_ID_Resp_good tr];
        disp(['Correct Resp at trial ' num2str(tr)])
        idx_Resp_g = [];
        idx_Resp_g = find(diff(evt_valve(idx_tr_st:idx_tr_end) & evt_lick(idx_tr_st:idx_tr_end))>0)+1; %Start resp & valve 
        idx_go = find(diff(evt_delay(idx_tr_st:idx_tr_end))<0) % End of delay 
        idx_puff = find(diff(evt_puff(idx_tr_st:idx_tr_end))>0)+1 % Start Puff
        
        if max(size(idx_Resp_g))>1
            idx_Resp_g = min(idx_Resp_g)
            
            disp(['Multiple Resp at trial ' num2str(tr)])
        end
        % idx within TRIAL
        idx_Resp_good = [idx_Resp_good idx_Resp_g]
        idx_Go_clean = [idx_Go_clean idx_go]% idx within TRIAL of the Response Lick that trigger rwd (the initiation)
        idx_Puff_clean = [idx_Puff_clean idx_puff]
    end
    disp(['nothing at trial ' num2str(tr)])
end
%%
disp(['Number of Clean Correct trials = ' num2str(max(size( tr_ID_Resp_good )))])
disp(['Percent of Clean Correct trials = ' num2str(100*(max(size( tr_ID_Resp_good ))/Ntrials)) '%'])

%% filter data
Nyquist= sr/2;
% [B,A] = butter(10,55/Nyquist, 'low'); % low pass filter (pass under 55Hz)
[B,A] = butter(3,[5/Nyquist 45/Nyquist]); % band pass filter (pass between 5-50Hz)
% [B,A] = butter(2,[55/Nyquist 100/Nyquist],'stop'); % band cut filter (cut between 55-100Hz)
dataf = filtfilt(B,A, data);

%% Select Data epoch of interest:
% ALL Correct trials

spectro_type = 'Puff_centered'

if spectro_type == 'Puff_centered'
    disp('Puff_centered')
    Nsec = 3 % second before and after the puff 
    dataf_Resp_good = [];
    idx_resp_g_trial =[];
    for tr_Rg =1:max(size(tr_ID_Resp_good)) % for each trial
        idx_resp_g_trial = idx_trial_start(tr_ID_Resp_good(tr_Rg)) + idx_Puff_clean(tr_Rg)
        dataf_Resp_good(tr_Rg,:)=dataf( idx_resp_g_trial - (Nsec*sr) : idx_resp_g_trial + (Nsec*sr));
        time_g_tr(tr_Rg,:)=time(idx_resp_g_trial - (Nsec*sr) : idx_resp_g_trial + (Nsec*sr));
        evt_puff_g_tr(tr_Rg,:)=evt_puff(idx_resp_g_trial - (Nsec*sr) : idx_resp_g_trial + (Nsec*sr));
        evt_delay_g_tr(tr_Rg,:)=evt_delay(idx_resp_g_trial - (Nsec*sr) : idx_resp_g_trial + (Nsec*sr));
        evt_lick_g_tr(tr_Rg,:)=evt_lick(idx_resp_g_trial - (Nsec*sr) : idx_resp_g_trial + (Nsec*sr));
    end
    
elseif spectro_type == 'GOcu_centered'
    disp('GOcue_centered')
    Nsec = 3 % second before and after the Go signal
    dataf_Resp_good = [];
    idx_resp_g_trial =[];
    for tr_Rg =1:max(size(tr_ID_Resp_good)) % for each trial
        idx_resp_g_trial = idx_trial_start(tr_ID_Resp_good(tr_Rg)) + idx_Go_clean(tr_Rg)
        dataf_Resp_good(tr_Rg,:)=dataf( idx_resp_g_trial - (Nsec*sr) : idx_resp_g_trial + (Nsec*sr));
        time_g_tr(tr_Rg,:)=time(idx_resp_g_trial - (Nsec*sr) : idx_resp_g_trial + (Nsec*sr));
        evt_puff_g_tr(tr_Rg,:)=evt_puff(idx_resp_g_trial - (Nsec*sr) : idx_resp_g_trial + (Nsec*sr));
        evt_delay_g_tr(tr_Rg,:)=evt_delay(idx_resp_g_trial - (Nsec*sr) : idx_resp_g_trial + (Nsec*sr));
        
    end
    
elseif spectro_type == 'Resp_centered'
    disp('RespL_centered')
    Nsec = 3 % second before and after the Resp
    dataf_Resp_good = [];
    idx_resp_g_trial =[];
    for tr_Rg =1:max(size(tr_ID_Resp_good)) % for each trial
        idx_resp_g_trial = idx_trial_start(tr_ID_Resp_good(tr_Rg)) + idx_Resp_good(tr_Rg)
        dataf_Resp_good(tr_Rg,:)=dataf( idx_resp_g_trial - (Nsec*sr) : idx_resp_g_trial + (Nsec*sr));
        %     time(tr_Rg,:)=idx_resp_g_trial - (Nsec*sr) : idx_resp_g_trial + (Nsec*sr);
    end
    
end
%% PLOT SPECTRO FOR EACH TRIALS
close all,
Pall=[];Sall=[];
for tr =1:max(size(tr_ID_Resp_good))
    figure,
    % Compute Spectrogram
    window=size(dataf_Resp_good,2)/(2*Nsec*6); noverlap=0; nfft= [10:1:80];
    [S,F,T,P] = spectrogram(dataf_Resp_good(tr,:),window,noverlap,nfft,sr,'yaxis');
    Pall(tr,:,:)= P;
    Sall(tr,:,:)= S;
    
    %     Plot raw Spectro for each trial
    %     figure, hold on
    %     subplot(8,1,4:6)
    %     colormap('jet')
    %     spectrogram(dataf_Resp_good(tr,:),window,noverlap,nfft,sr,'yaxis');
    %     caxis([-0 50])
    %     title(['Correct trial #' num2str(tr_ID_Resp_good(tr)) ])
    
    %     Plot Normalized Spectro for each trial
    nP=[];
    for i=1 : max(size(T))
        totP = sum(P(:,i));
        nP(:,i) = P(:,i) ./ totP; 
    end
    nPall(tr,:,:)= nP;
    
    
    % nP = normalized P 
    subplot(10,1,1:4),  hold on,
    surf(T,F,nP,'EdgeColor','none');
    axis xy; axis tight;
    colormap(jet);
    colorbar;
    view(0,90);
    ylim([10 60])
    %     xlabel('Time (sec)')
    ylabel('Frequencies (Hz)')
    title (['Normalized Spectro for trial #' num2str(tr_ID_Resp_good(tr)) ])
    
    
    subplot(10,1,5:6), hold on,
    plot(time_g_tr(tr,:)-time_g_tr(tr,1), dataf_Resp_good(tr,:))
    plot(time_g_tr(tr,:)-time_g_tr(tr,1), evt_puff_g_tr(tr,:)*100,'r' )
    plot(time_g_tr(tr,:)-time_g_tr(tr,1), evt_delay_g_tr(tr,:)*100,'g')
    xlim([0 6])
    ylim([-200 200])
    colorbar;
    
    % P = not normalized 
    subplot(10,1,8:10),  hold on
    z = log(P)*5; %e.g. 37x17
    t = T;      %e.g  1x17
    x = F;      %e.g. 37x1
    surf(t,x,z,'EdgeColor','none');
    axis xy; axis tight;
    colormap(jet);
    colorbar
    caxis([-10 10])
    view(0,90);
    xlabel('Time (sec)')
    
end
%% AVERAGED Normalized Spectrogram (Normalized to the total power at each time bin)
figure,

nP=[];
nPmean= squeeze(mean(nPall,1));
for i=1 : max(size(T))
    totP = sum(nPmean(:,i));
    nP(:,i) = nPmean(:,i) ./ totP;
end

subplot(10,1,1:4),  hold on,
surf(T,F,nP,'EdgeColor','none');
axis xy; axis tight;
colormap(jet);
colorbar;
view(0,90);
ylim([10 60])
%     xlabel('Time (sec)')
ylabel('Frequencies (Hz)')
title (['Averaged Normalized Spectro'])

evt_ALL_licks = evt_puff_g_tr(1,:)*0;
subplot(10,1,5:6), hold on,
for tr= 1:max(size(tr_ID_Resp_good))
    evt_ALL_licks = evt_ALL_licks + evt_lick_g_tr(tr,:);
end


plot(time_g_tr(tr,:)-time_g_tr(tr,1), evt_ALL_licks*100,'k');
plot(time_g_tr(tr,:)-time_g_tr(tr,1), evt_puff_g_tr(tr,:)*500,'r' );
plot(time_g_tr(tr,:)-time_g_tr(tr,1), evt_delay_g_tr(tr,:)*500,'g');
xlim([0 6]);
%     ylim([-200 200]);
legend('Licks','Puff','delay','Location','NorthWest')
c=colorbar;
c.Visible = 'off'


subplot(10,1,8:10),  hold on
Pmean= squeeze(mean(Pall,1));
z = log(Pmean)*5; %e.g. 37x17;
t = T;      %e.g  1x17
x = F;      %e.g. 37x1
surf(t,x,z,'EdgeColor','none');
axis xy; axis tight;
colormap(jet);
colorbar
caxis([-10 10])
view(0,90);
xlabel('Time (sec)')


% plot raw traces (filtered data)

% figure, hold on
% for i=1:10;
%     plot(dataf_Resp_good(i,:))
% end
% legend

% %% replicate plot with surf
% figure,
% z = log(P)*5; %e.g. 37x17
% t = T;      %e.g  1x17
% x = F;      %e.g. 37x1
% surf(t,x,z,'EdgeColor','none');
% axis xy; axis tight;
% colormap(jet);
% colorbar
% caxis([0 50])
% view(0,90);

% %% plot average
% Pmean= squeeze(mean(Pall,1))
%
% figure,
% z = log(abs(Pmean))*5;
% t = T; %1x17
% x = F; %37x1
% surf(t,x,z,'EdgeColor','none');
% axis xy; axis tight;
% colormap(jet);
% colorbar
% caxis([0 50])
% view(0,90);
% xlabel('Time (sec)')
% ylabel('Frequencies (Hz)')
% title (['Averaged Spectrogram Centered on the Resp (+/-' num2str(Nsec) 'sec)'])
% xlim([1.4 4.6])
% ylim([10 50])



