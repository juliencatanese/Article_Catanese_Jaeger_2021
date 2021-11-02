% plot_PSDsub_vs_PSDraw_1Chan_JCscript
% Plot PSDsub vs PSDraw for 1 Chanel manually selected.
% Also plot Normalized PSD (PSDnorm = PSD x 1/freq^2).
% Using pwelch (Matlab) , Hanning Windows, wsize = 2^13;
% Written by Catanese J., 13 August 2017, Jaegerlab, Emory.
% Last updated on 11/12/2018

close all; clear all,

% colormap summer
%% Define your Channel of interest
Sh=2;
Ch=4;
ChanID= ['S' num2str(Sh) 'Ch' num2str(Ch)];
Ref1ID= ['S' num2str(Sh) 'Ch' num2str(Ch+1)];
Ref2ID= ['S' num2str(Sh) 'Ch' num2str(Ch-1)];
load('info.mat');
MouseID=info.info_notes.MouseID
Day=info.info_notes.Day

%% Get DATA
load([ChanID '_raw.mat']); data_Ori = data; clear data;
load([Ref1ID '_raw.mat']); data_Ref01 = data; clear data;
load([Ref2ID '_raw.mat']); data_Ref02 = data; clear data;

%% restrict data to Eoch of interest
load('Epochs_pre_post_task_st_end.mat','idx_preStim_st_end', 'idx_task_st_end', 'idx_postStim_st_end', 'idx_postRest_st_end', 'idx_preRest_st_end')
figure(1), close(1),
for Epo=1:3
    if Epo==1
        EpochID = 'preRest'
        idx_St_End= idx_preRest_st_end;
    elseif Epo==2
        EpochID= 'task'
        idx_St_End= idx_task_st_end;
    elseif Epo==3
        EpochID= 'postRest'
        idx_St_End= idx_postRest_st_end;
    end
    data_raw = data_Ori(idx_St_End(1):idx_St_End(2));
    data_Ref1 = data_Ref01(idx_St_End(1):idx_St_End(2));
    data_Ref2 = data_Ref02(idx_St_End(1):idx_St_End(2));
    data_sub = data_raw - mean([data_Ref1; data_Ref2]);
    
    %% Compute PSD
    wsize = 2^13      % (Hanning windows size to compute pwelch)
    
    for raw_sub=1:2
        if raw_sub==1
            datID='raw'
            [PSD,F] = pwelch(data_raw,hanning(wsize),wsize/2,2*wsize,sr); % sr = sampling rate
        else
            datID='sub'
            [PSD,F] = pwelch(data_sub,hanning(wsize),wsize/2,2*wsize,sr); % sr = sampling rate
        end
        
        %% Normalization
        PSDnorm=PSD;
        for i=1:max(size(PSDnorm))
            PSDnorm(i,:)= PSDnorm(i,:) .* F(i)^2;
        end
        
        %% Plot Power Spectrum (per shanks)
        if raw_sub==2; col=[0.3 0.3 0.3]; else ; col=[0.7 0.7 0.7]; end

        figure(1),
        subplot(2,3,Epo), hold on,
        hold on, plot(F,10*log10(PSDnorm),'Color', col ,'LineWidth',1.5)
        gcf, xlim([0 100]);
        gcf, ylim([10 80]);
        gcf, xlabel('Freq (Hz)')
        gcf ,ylabel('Power (dB)')
        gcf, title([ EpochID ])
        
        
        subplot(2,3,Epo+3), hold on,
        plot(F,10*log10(PSDnorm),'Color', col ,'LineWidth',1.5)
        gcf, xlim([0 40]);
        gcf, ylim([10 80]);

        gcf, xlabel('Freq (Hz)')
        gcf ,ylabel('Power (dB)')
        gcf, title([MouseID ' ' Day ' PSDn'])
        
    end
end
% legend
legend([ChanID 'raw'] , [ChanID 'sub'] ,'Location','Best')
% SAVING
mkdir('.\LFP\')
saveas(gcf, ['.\LFP\' ChanID '_PSDn_Epochs_pre_task_post'],'png')
saveas(gcf, ['D:\JC_Figures\LFP\PSD\' MouseID '_' Day  '_' ChanID '_PSDn__Epochs_pre_task_post'],'png')