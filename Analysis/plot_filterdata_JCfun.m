function plot_filterdata_JCfun(Shank, Chanel, minute_start, minute_end)
% function plot_rawdata_JCfun(SXChX, minute_start, minute_end)
% define which minute interval of raw data you want to visualized (from 1 to 20min or more see maxtime)
% JC last updated 9/24/2018
load(['S' num2str(Shank) 'Ch' num2str(Chanel) '_raw.mat']); % data in microVolt 
% minute_start = 3 
% minute_end = 4

% load('S1Ch7.mat')
load('info.mat');
load('time.mat');

% reduce plot to a certain time interval
sr=info.info_freq_parameters.amplifier_sample_rate;
maxtime = max(time)/60; % in minute
disp(['max time is ' num2str(maxtime) 'minutes'])

%% Filter data for spikes
[B,A] = butter(2,250/10000,'high'); % high pass filter
data_filt = filtfilt(B,A, data);

%% threshold
thr = -abs(3*std(data_filt)); % threshold is -3std from filtered signal.
thr_bool = data_filt < thr; % valeur de voltage inferieur au threshold
thr_line = thr*ones(1,max(size(time))); % trh line for figure

%% PLOT 
% figure, 
plot(time(sr*60*minute_start:sr*60*minute_end), data_filt(sr*60*minute_start:sr*60*minute_end))
title(['S' num2str(Shank) 'Ch' num2str(Chanel)])
ylabel('microV'); 
xlabel('time (s)');
