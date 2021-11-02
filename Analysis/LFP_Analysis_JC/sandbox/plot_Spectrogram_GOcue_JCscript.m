% 



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
idx_go =   find(diff(evt_delay)<0); %idx within the session of the end of the delay (=Go Signal)

evt_valve = evt_rwd_L + evt_rwd_R;  
evt_lick = evt_lick_L + evt_lick_R; 
% idx_Resp_good = find(diff(evt_valve & evt_lick)>0)+1; % idx within sessions of the start of the Response Lick that trigger rwd   

%% Trial selection : 
% correct trials
tr_ID_Resp_good = []; 
idx_Resp_good= []; 
for tr=1:75-1
    % Collect data for each trial
   if sum(evt_valve(idx_trial_start(tr):idx_trial_end(tr)) & evt_lick(idx_trial_start(tr):idx_trial_end(tr)))>0
    tr_ID_Resp_good = [tr_ID_Resp_good tr];
    disp(['Correct Resp at trial ' num2str(tr)])
    idx_Resp_g = []; 
    idx_Resp_g = find(diff(evt_valve(idx_trial_start(tr):idx_trial_end(tr)) & evt_lick(idx_trial_start(tr):idx_trial_end(tr)))>0)+1; 
        if max(size(idx_Resp_g))>1 
            idx_Resp_g = min(idx_Resp_g)
            disp(['Multiple Resp at trial ' num2str(tr)])
        end
    idx_Resp_good = [idx_Resp_good idx_Resp_g]    
   end
            disp(['nothing at trial ' num2str(tr)])

end

%% filter data
Nyquist= sr/2;
% [B,A] = butter(2,[5/Nyquist 200/Nyquist]); % low pass filter 150Hz
[B,A] = butter(2,[50/Nyquist 180/Nyquist],'stop'); % low pass filter 150Hz

dataf = filtfilt(B,A, data);

%% Spectro around GO for each trials
Nsec = 5 % second before and after the Go signal
dataf_go = [];
for tr=1:75
    % Collect data for each trial
    dataf_go(tr,:)=dataf(idx_go(tr)-(Nsec*sr):idx_go(tr)+(Nsec*sr));
end
%%
Pall=[];Sall=[];
for tr=1:75
    % Compute Spectrogram
    window=size(dataf_go,2)/(2*Nsec*5); noverlap=round(window/4); nfft= [15:1:80];
    [S,F,T,P] = spectrogram(dataf_go(tr,:),window,noverlap,nfft,sr,'yaxis');
    Pall(tr,:,:)= P;
    Sall(tr,:,:)= S;
    
end
%% plot example
close all, 

figure,
colormap('jet')
% tr=49% indicate the trial you want to see: e.g. trial #1
spectrogram(dataf_go(tr,:),window,noverlap,nfft,sr,'yaxis');
caxis([-0 30])

%% replicate plot with surf
% z = randn(100,100);
% t = 1:100;
% x = 1:100;
% surf(t,x,abs(z),'EdgeColor','none');
% axis xy; axis tight; colormap(jet); view(0,90);
figure, 
z = log(P)*5;
t = T; %1x17
x = F; %37x1
surf(t,x,z,'EdgeColor','none');
axis xy; axis tight; 
colormap(jet); 
colorbar
caxis([0 60])
view(0,90);


%% plot average

Pmean= squeeze(mean(Pall,1))

figure, 
z = log(abs(Pmean))*5;
t = T; %1x17
x = F; %37x1
surf(t,x,z,'EdgeColor','none');
axis xy; axis tight; 
colormap(jet); 
colorbar
caxis([-10 50])
view(0,90);
xlabel('Time (sec) [center = 5s = GO]')
ylabel('Frequencies (Hz)')
title ('Averaged Spectrogram Centered on the GO sound (+/- 5sec)')



