function plot_rawdata_JCfun(Shank, Chanel, minute_start, minute_end)
% function plot_rawdata_JCfun(SXChX, minute_start, minute_end)
% example of input arg values: Shank= 1;  Chanel=1; , minute_start=10; , minute_end=11 ; 
% define which minute interval of raw data you want to visualized (from 1 to 20min or more see maxtime)
% Written by Julien Catanese 
% last update JC 9/26/2018

% load(['S' num2str(Shank) 'Ch' num2str(Chanel) '_raw.mat']);
load(['S' num2str(Shank) 'Ch' num2str(Chanel) '_sub.mat']);

% minute_start = 3 
% minute_end = 4

% load('S1Ch7.mat')
load('info.mat');
load('time.mat');

SXChX = [];
SXChX = data; 

SFreq=info.info_freq_parameters.amplifier_sample_rate
sr=SFreq; 

maxtime = max(time)/60; % in minute
disp(['max time is ' num2str(maxtime) 'minutes'])

% figure, 
plot(time(sr*60*minute_start:sr*60*minute_end),SXChX(sr*60*minute_start:sr*60*minute_end))
title(['S' num2str(Shank) 'Ch' num2str(Chanel)])
ylabel('microV'); 
xlabel('time (s)');
