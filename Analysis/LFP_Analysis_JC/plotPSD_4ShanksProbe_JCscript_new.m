% Script plotPSD_4ShanksProbe_fun_JC(dataFolder, figFolder, wsize, Xlim)
% SubPlot(2,4) the PSD of each channel (color legend) per shank (4 columns)
% the 2nd row are Normalized PSD (PSDnorm = PSD x 1/freq^2).
% wsize = 2^13      % (Hanning windows size to compute pwelch)
% Xlim = 40  %(figure Xaxis ploted from 0 Hz to 40 Hz)
% Written by Catanese J., 13 August 2017, Jaegerlab, Emory.
% Last updated on 11/12/2018

close all; clear all,
colormap summer

%% Define Mouse and Chanels of Interest
MouseID = 'JCVGAT15'
ChanID={'S1Ch6','S2Ch1', 'S3Ch4', 'S4Ch1'}
tag_ChanID_all = []

%% Loop For each Channel
for nCh= 1:max(size(ChanID))
    %  load data
    datID='sub'
    load([ChanID{nCh} '_sub.mat'])
    disp(ChanID{nCh})
    
    %% Compute PSD
    wsize = 2^13      % (Hanning windows size to compute pwelch)
    [PSD,F] = pwelch(data,hanning(wsize),wsize/2,2*wsize,sr); % sr = sampling rate
    %% Normalization
    PSDnorm=PSD;
    for i=1:max(size(PSDnorm))
        PSDnorm(i,:)= PSDnorm(i,:) .* F(i)^2;
    end
    
    %% plot
    hold on,
    plot(F,10*log10(PSDnorm),'LineWidth',1.5)  %hold on;\
    gcf, xlim([0 100]);
    gcf, xlabel('Freq (Hz)')
    gcf ,ylabel('Power (dB)')
    gcf, title(['PSD normalized'])
    gcf, legend(ChanID,'Location','southwest')
    tag_ChanID_all = [tag_ChanID_all  ChanID{nCH}]
end
%% SAVING
saveas(gcf, ['.\PSD_multiChan_' tag_ChanID_all ],'png')
saveas(gcf, ['D:\JC_Figures\LFP\PSD\PSD_multiChan_' tag_ChanID_all],'png')
