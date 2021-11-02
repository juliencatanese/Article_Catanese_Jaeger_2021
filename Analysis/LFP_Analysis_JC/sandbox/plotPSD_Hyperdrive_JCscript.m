% Script plotPSD_Hyperdrive_JCscript
% SubPlot(2,4) the PSD of each channel (color legend) per shank (4 columns)
% the 2nd row are Normalized PSD (PSDnorm = PSD x 1/freq^2).
% wsize = 2^13      % (Hanning windows size to compute pwelch)
% Xlim = 40  %(figure Xaxis ploted from 0 Hz to 40 Hz)
% by Catanese J., 13 August 2017, Jaegerlab, Emory.

close all; clear all,

%% arguments:
dataFolder =  'D:\JC_Data\Chronic\HDr_JC02\Day1_Nochoicetask_500msDelay\'    % 'D:\JC_Data\Acute\Awake\JCTITL01\Insertion1\'
figFolder =  'D:\JC_Figures\Chronic\HDr_JC02\Day1_Nochoicetask_500msDelay\'  % 'D:\JC_Figures\Acute\Awake\JCTITL01\Insertion1\'
mkdir(figFolder)

% wsize = 2^13      % (Hanning windows size to compute pwelch)
% Xlim = 40  %(figure Xaxis ploted from 0 Hz to 40 Hz)

d=dir([dataFolder '*.rhd'])

%% Loop through Files in Folder
for nfile = 1:size(d,1)
    nfile
    File = d(nfile).name
    Folder = dataFolder
    
    %% Get DATA
    read_Intan_RHD2000_file_JC(Folder, File) % type 'whos' to see the resulting variables
    
    %% Compute PSD
    
    Fs = 20000;
    wsize = 2^12      % (Hanning windows size to compute pwelch)
    [PSD,F] = pwelch(amplifier_data',hanning(wsize),wsize/2,2*wsize,Fs);
    %% Normalization
    PSDnorm=PSD;
    for i=1:max(size(PSDnorm))
        PSDnorm(i,:)= PSDnorm(i,:) .* F(i)^2;
    end
    
    %% Plot Power Spectrum (per All Channels)
    
    Xlim = 60  %(figure Xaxis ploted from 0 Hz to 40 Hz)
    
    
    figure(nfile),
    subplot(2,2,1),
    plot(F,10*log10(PSD(:,:)),'LineWidth',1.5); %hold on;
    set(gca,'YLim',[10 90], 'XLim',[0 300],'FontSize',12); grid on;
%     legend({'ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'ch7', 'ch8'},'Location','Northeast')
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    gcf, title([File(1:7) ' PSD ']);
    
    subplot(2,2,2),
    plot(F,10*log10(PSD(:,:)),'LineWidth',1.5); %hold on;
    set(gca,'YLim',[10 90], 'XLim',[0 Xlim],'XTick',0:Xlim/4:Xlim,'FontSize',12); grid on;
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    gcf, title([File(1:7) ' PSD ']);
    
    
    subplot(2,2,3),
    plot(F,10*log10(PSDnorm(:,:)),'LineWidth',1.5); %hold on;
    set(gca,'YLim',[50 90], 'XLim',[0 300],'FontSize',12); grid on;
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    gcf, title([File(1:7) ' PSD Norm'])
    
    subplot(2,2,4),
    plot(F,10*log10(PSDnorm(:,:)),'LineWidth',1.5); %hold on;
    set(gca,'YLim',[50 90], 'XLim',[0 Xlim],'XTick',0:Xlim/4:Xlim,'FontSize',12); grid on;
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    gcf, title([File(1:7) ' PSD Norm'])
    
    
    
    %% SAVING
%     saveas(gcf, [ figFolder 'PSD-' File(1:7) '#' num2str(nfile)  '-hanWindw' num2str(wsize) '-Xlim' num2str(Xlim)  'Hz'   '.m'],'m')
    saveas(gcf, [ figFolder 'PSD-' File(1:7) '#' num2str(nfile) '-hanWindw' num2str(wsize) '-Xlim' num2str(Xlim)  'Hz'  '.png'],'png')
end


