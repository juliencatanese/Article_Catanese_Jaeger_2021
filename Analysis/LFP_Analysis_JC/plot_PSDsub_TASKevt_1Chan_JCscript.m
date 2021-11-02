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
idx_St_End= idx_task_st_end;
data_raw = data_Ori(idx_St_End(1):idx_St_End(2));
%% Substract local Reference (mean of 2 closest channels)
data_Ref1 = data_Ref01(idx_St_End(1):idx_St_End(2));
data_Ref2 = data_Ref02(idx_St_End(1):idx_St_End(2));
data_sub = data_raw - mean([data_Ref1; data_Ref2]);
figure(1), close(1),
%%
load('info.mat');
load('time.mat');
load('Ntrial_type.mat', 'trial')

for Epo=1:4
    
    if Epo==1
        EpochID = 'Sample'
        load('evt.mat','evt_puff_R', 'evt_puff_L');
        evt_puff = evt_puff_L + evt_puff_R; clear evt_puff_R evt_puff_L;
        evt= evt_puff(idx_St_End(1):idx_St_End(2));
        evt_st = find(diff(evt)>0); % the start of the puff
        evt_end = find(diff(evt)<0);
        
    elseif Epo==2
        EpochID= 'Delay'
        load('evt.mat','evt_delay');
        evt= evt_delay(idx_St_End(1):idx_St_End(2));
        evt_st = find(diff(evt)>0); % the start of the delay
        evt_end = find(diff(evt)<0);
        
    elseif Epo==3
        EpochID= 'Response' 
        evt= evt_delay(idx_St_End(1):idx_St_End(2));
        evt_st = find(diff(evt)<0); % the end of the delay
        evt_end = find(diff(evt)<0)+0.750*sr;

    elseif Epo==4
        EpochID= 'Pre Lick'
        load('evt.mat','evt_trial','evt_lick_L', 'evt_lick_R');      
        evt_lick= evt_lick_L + evt_lick_R; clear evt_lick_R evt_lick_L;
        evt_lick = evt_lick(idx_St_End(1):idx_St_End(2));
        evt_tr= evt_trial(idx_St_End(1):idx_St_End(2));
        tr_st = find(diff(evt)>0);  % the start of the delay
        tr_end = find(diff(evt_tr)>0); % the end of the trial
        tr_end=tr_end(2:end);
        firstLick = [];
        for tr=1:(trial.Ntrial)
%               diff([tr_st(tr) tr_end(tr)])
            A = min(find(evt_lick(tr_st(tr):tr_end(tr))))
            if isempty(A)
                idxFlick = 0;
            else
                idxFlick= -1 + tr_st(tr) + min(find(evt_lick(tr_st(tr):tr_end(tr))))    ;
            end
            firstLick = [firstLick  idxFlick];
            evt_st = firstLick-0.750*sr  ;
            evt_end = firstLick        ;
        end
    end 
        %% Separates trials types
        idx_corr = sort([trial.idx_correct_L trial.idx_correct_R]);
        
        idx_evt_st=evt_st(idx_corr);
        idx_evt_end=evt_end(idx_corr);
        
        
        %% Compute PSD restricted to a given Task Epoch (eg delay)
        
        for raw_sub=1:2
            if raw_sub==1
                datID='raw'
                data = [];
                for tr=1:max(size(idx_corr))
                    data = [data data_raw(idx_evt_st(tr):idx_evt_end(tr))];
                end
            else
                datID='sub'
                data = [];
                for tr=1:max(size(idx_corr))
                    data = [data data_sub(idx_evt_st(tr):idx_evt_end(tr))];
                end
                
            end
            
            wsize = 2^13      % (Hanning windows size to compute pwelch)
            [PSD,F] = pwelch(data,hanning(wsize),wsize/2,2*wsize,sr); % sr = sampling rate
            
            %% Normalization
            PSDnorm=PSD;
            for i=1:max(size(PSDnorm))
                PSDnorm(i,:)= PSDnorm(i,:) .* F(i)^2;
            end
            
            %% Plot Power Spectrum (per shanks)
            if raw_sub==2; col=[0.3 0.3 0.3]; else ; col=[0.7 0.7 0.7]; end
            
            figure(1),
            subplot(2,4,Epo), hold on,
            hold on, plot(F,10*log10(PSDnorm),'Color', col ,'LineWidth',1.5)
            gcf, xlim([0 100]);
            gcf, ylim([10 80]);
            gcf, xlabel('Freq (Hz)')
            gcf ,ylabel('Power (dB)')
            gcf, title([ EpochID ])
            
            
            subplot(2,4,Epo+4), hold on,
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
    saveas(gcf, ['.\LFP\' ChanID '_PSDn_TASK_evts'],'png')
    saveas(gcf, ['D:\JC_Figures\LFP\PSD\' MouseID '_' Day  '_' ChanID '_PSDn_TASK_evts'],'png')