% plot_PSDsub_vs_PSDraw_1Chan_JCscript
% Plot PSDsub vs PSDraw for 1 Chanel manually selected.
% Also plot Normalized PSD (PSDnorm = PSD x 1/freq^2).
% Using pwelch (Matlab) , Hanning Windows, wsize = 2^13;
% Written by Catanese J., 13 August 2017, Jaegerlab, Emory.
% Last updated on 11/12/2018

close all; clear all,
colormap summer
%%
load('info.mat')
MouseID=info.info_notes.MouseID
Day=info.info_notes.Day
ChanID='S1Ch6'

%% Get DATA

for dat=1:2;
    if dat==1
        datID='sub';
        load([ChanID '_sub.mat'])
    else
        datID='raw';
        load([ChanID '_raw.mat'])
    end
    
    disp(ChanID)
    %% Compute PSD
    wsize = 2^13      % (Hanning windows size to compute pwelch)
    [PSD,F] = pwelch(data,hanning(wsize),wsize/2,2*wsize,sr); % sr = sampling rate
    %% Normalization
    PSDnorm=PSD;
    for i=1:max(size(PSDnorm))
        PSDnorm(i,:)= PSDnorm(i,:) .* F(i)^2;
    end
    
    %% Plot Power Spectrum (per shanks)
    
    %         close all,
    figure(1)
    if dat==1; col=[0.3 0.3 0.3]; else ; col=[0.7 0.7 0.7]; end
    gcf, hold on,
    subplot(2,2,1), hold on,
    plot(F,10*log10(PSDnorm),'Color', col ,'LineWidth',1.5)
    gcf, xlim([0 100]);
    gcf, xlabel('Freq (Hz)')
    gcf ,ylabel('Power (dB)')
    gcf, title(['PSD normalized ' ChanID ' ' MouseID])
    legend('PSDnorm = PSD.F^2)','Location','southwest')
    
    subplot(2,2,2), hold on,
    plot(F,10*log10(PSD),'Color', col ,'LineWidth',1)
    gcf, xlim([0 100]);
    gcf, xlabel('Freq (Hz)')
    gcf ,ylabel('Power (dB)')
    gcf, title(['PSD ' ChanID ' ' MouseID])
    if dat==1
        legend(['PSD ' datID])
    elseif dat==2
        legend('PSD sub', ['PSD' datID])
    end
    
    subplot(2,2,3), hold on,
    plot(F,10*log10(PSDnorm),'Color', col ,'LineWidth',1.5)
    gcf, xlim([0 40]);
    gcf, xlabel('Freq (Hz)')
    gcf ,ylabel('Power (dB)')
    
    subplot(2,2,4), hold on,
    plot(F,10*log10(PSD),'Color', col ,'LineWidth',1)
    gcf, xlim([0 40]);
    gcf, xlabel('Freq (Hz)')
    gcf ,ylabel('Power (dB)')
    mkdir('.\LFP\')
    saveas(gcf, ['.\LFP\' ChanID '_PSDn_PSD'],'png')
    saveas(gcf, ['D:\JC_Figures\LFP\PSD\' MouseID '_' Day  '_' ChanID '_PSDn_PSD_sub_vs_raw'],'png')
end