% Script plot_PSD_MultiChan_JCscript
% Plot and compare normPSD of each channel in the list (e.g. ChanID={'S1Ch6','S2Ch1', 'S3Ch4', 'S4Ch1'})
% Written by Catanese J., 13 August 2017, Jaegerlab, Emory.
% Last updated on 11/12/2018

close all; clear all;
colormap summer;

load('info.mat');
MouseID = info.info_notes.MouseID

%% Define Chanels List
ChanID={'S1Ch1', 'S1Ch2', 'S1Ch3','S1Ch5', 'S1Ch6', 'S1Ch7', 'S1Ch8'}
% ChanID={'S3Ch1','S3Ch4', 'S3Ch8'}
% ChanID={'S4Ch1', 'S4Ch4','S4Ch5', 'S4Ch6', 'S4Ch7', 'S4Ch8'}
% % ChanID={'S1Ch1','S2Ch1', 'S3Ch1', 'S4Ch1'; 'S1Ch8','S2Ch8', 'S3Ch8', 'S4Ch7'}
% ChanID={'S1Ch1', 'S1Ch3','S1Ch5', 'S1Ch7'; 'S2Ch1', 'S2Ch3', 'S2Ch5', 'S2Ch7'}
datID='sub'
tag_ChanID_all = [];

%% Loop For each Channel
for nRow= 1:min(size(ChanID))
    subplot(min(size(ChanID)),1,nRow);
    
    for nCh= 1:max(size(ChanID))
        %  load data
        
        load([ChanID{nRow,nCh} '_' datID '.mat']);
        disp(ChanID{nRow,nCh});
        
        %% Compute PSD
        wsize = 2^13 ;     % (Hanning windows size to compute pwelch)
        [PSD,F] = pwelch(data,hanning(wsize),wsize/2,2*wsize,sr); % sr = sampling rate
        %% Normalization
        PSDnorm=PSD;
        for i=1:max(size(PSDnorm));
            PSDnorm(i,:)= PSDnorm(i,:) .* F(i)^2;
        end
        
        %% plot
        hold on,
        plot(F,10*log10(PSDnorm),'LineWidth',1.5);  %hold on;\
        gcf, xlim([0 70]);
        gcf, xlabel('Freq (Hz)');
        gcf ,ylabel('Power (dB)');
        gcf, title(['PSD norm ' MouseID]);
        gcf, legend(ChanID{nRow,:},'Location','southwest');
        tag_ChanID_all = [tag_ChanID_all  ChanID{nRow,nCh}];
        
    end
end
%% SAVING
saveas(gcf, ['.\PSD_multiChan_' tag_ChanID_all ],'png');
saveas(gcf, ['D:\JC_Figures\LFP\PSD\PSD_multiChan_' MouseID '_' tag_ChanID_all],'png');
