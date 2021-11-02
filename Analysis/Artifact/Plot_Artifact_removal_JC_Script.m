% Plot_Artifact_removal_JC_Script

% S2Ch7 = ref (the one with the less number of Spikes after AutoClustering) 
clear all, 

min_st= 5; 
min_end=min_st+0.1; 
RefChan='S2Ch7'
Sh=1; 
Ch=2; 
load([RefChan '_raw.mat']); ref=data(sr*60*min_st:sr*60*min_end); clear data; 
load(['S' num2str(Sh) 'Ch' num2str(Ch) '_raw.mat']); ch1=data(sr*60*min_st:sr*60*min_end); clear data;
[B,A] = butter(2,600/10000,'high'); % high pass filter

filt_ref = filtfilt(B,A, ref);
filt_ch1 = filtfilt(B,A, ch1);

load('time.mat'); t=time(sr*60*min_st:sr*60*min_end);

%% Compare plot 
close all, 

figure,

subplot(3,1,1), hold on, 
% plot raw data  
plot(t, ch1),
plot(t, ref)
legend('ch1', 'ref', 'Location','best' )
ylim([-500 500])

subplot(3,1,2), hold on, 
% substraction on Raw data.   
ch1_rawsub= ch1 - ref; 
plot(t, ch1)
plot(t, ch1_rawsub)
legend('ch1_raw', 'ch1_raw - ref_raw', 'Location','best' )  
ylim([-500 500])

subplot(3,1,3), hold on, 

% filter without substraction
plot(t, filt_ch1)
% filtering before substraction  
ch1_filtsub = filt_ch1 - filt_ref; 
plot(t, ch1_filtsub)  
% filtering after substraction 
ch1_rawsub_filt = filtfilt(B,A, ch1_rawsub);
plot(t, ch1_rawsub_filt+300) 

ylim([-500 800])
legend( ' filter without substraction', 'filtering before substraction  ', 'filtering after substraction ' , 'Location','best' )



%% threshold
thr = -abs(2*std(ch1_rawsub_filt)); % threshold is -3std from filtered signal.
thr_line1 = thr*ones(1,max(size(t))); % trh line for figure

thr = -abs(20*std(ch1_rawsub_filt)); % threshold is -3std from filtered signal.
thr_line2 = thr*ones(1,max(size(t))); % trh line for figure

% thr_bool = ch1_rawsub_filt < thr; % valeur de voltage inferieur au threshold

figure, hold on, 
plot(t, ch1_rawsub_filt) 
plot(t, thr_line1)
plot(t, thr_line2)
ylim([-800 500])

