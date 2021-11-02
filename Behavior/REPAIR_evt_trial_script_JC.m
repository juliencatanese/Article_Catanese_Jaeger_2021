%  behavior_master_script_JC
%  play  this for each folder contening "board-DIN-00.dat"
%  "board-DIN-01.dat" ...
%  JCatanese 2017 in  JaegerLab

%%
clear all, close all,
%% To run a loop over folder
MouseID = 'JCVGAT1';
FolderID = dir(['D:\JC_Analysis\' MouseID '*\*1807*']);
% FolderID = dir(['D:\JC_Data\Acute\Awake\' MouseID '\w*d*S3Lv*1806*']);
% FolderID =dir(['D:\JC_Data\Acute\Awake\' MouseID '\w*d*rec*notri*\z*'])
%%

for nf=1:max(size(FolderID))
    cd([FolderID(nf).folder '\' FolderID(nf).name])
    disp([num2str(nf) '/' num2str(max(size(FolderID))) '   MouseID =' MouseID  '   FolderID = ' FolderID(nf).name])
    
    % convert .dat into .mat
    if isempty(dir('*.mat'))
        dat2mat_JC_Script
    else
        disp('already computed')
    end
    
    %% load evts.mat
    load('info.mat');
    load('evt.mat');
    load('time.mat');
    sr=info.info_freq_parameters.board_dig_in_sample_rate;
    evt_puff = evt_puff_L + evt_puff_R;
    maxtime=time(end);
    midtime=time(end/2);
    
    %% plot a figure to represent the type of trials
    if isempty(dir('evtfig.png'))
        figure, hold on
        subplot(2,1,1), hold on,
        plot(time(((maxtime*(1/4))*sr:((maxtime*(1/4))+50)*sr)), evt_trial(((maxtime*(1/4))*sr:((maxtime*(1/4))+50)*sr)),'b')
        plot(time(((maxtime*(1/4))*sr:((maxtime*(1/4))+50)*sr)), evt_puff(((maxtime*(1/4))*sr:((maxtime*(1/4))+50)*sr)).*0.5,'r')
        plot(time(((maxtime*(1/4))*sr:((maxtime*(1/4))+50)*sr)), evt_delay(((maxtime*(1/4))*sr:((maxtime*(1/4))+50)*sr)).*0.3,'g')
        legend('trigtrial','puffL','delay')
        xlabel('1min time (s)')
        subplot(2,1,2), hold on,
        plot(time(((maxtime*(3/4))*sr:((maxtime*(3/4))+50)*sr)), evt_trial(((maxtime*(3/4))*sr:((maxtime*(3/4))+50)*sr)),'b')
        plot(time(((maxtime*(3/4))*sr:((maxtime*(3/4))+50)*sr)), evt_puff(((maxtime*(3/4))*sr:((maxtime*(3/4))+50)*sr)).*0.5,'r')
        plot(time(((maxtime*(3/4))*sr:((maxtime*(3/4))+50)*sr)), evt_delay(((maxtime*(3/4))*sr:((maxtime*(3/4))+50)*sr)).*0.3,'g')
        title('before correction')
        saveas(gcf,'evtfig','png')
    else
        disp('already ploted evtfig.png')
    end
    
    
    %%  re-create the evt_trial.mat
    idx_start_puff = find(diff(evt_puff)>0)+1;
    idx_end_puff = find(diff(evt_puff)<0)+1;
    idx_end_trigtrial= idx_start_puff-2;
    idx_start_trigtrial= idx_end_trigtrial-(sr*0.01); % 10ms TTL for trigger trial just before the puff.
    
    
    evt_trial(:,:)=0;
    for Nsr=1:sr*0.01 %10ms
        evt_trial(idx_start_trigtrial+Nsr) = 1;
        % evt_trial(idx_start_trigtrial) = 1;
    end
    
    for tr=1:max(size(idx_start_trigtrial))
        if evt_delay(idx_end_puff(tr)+0.05*sr) == 0
            disp(['ATTENTION: Incompplete trial #' num2str(tr)])
            for Nsr=1:sr*0.01 %10ms
                evt_delay(idx_end_puff(tr)+Nsr+1) = 1;
            end
        end
    end
    
    save('evt.mat','evt_trial', 'evt_opto', 'evt_delay', 'evt_puff_L', 'evt_puff_R', 'evt_rwd_L', 'evt_rwd_R', 'evt_lick_L', 'evt_lick_R')
    
    %% test to plot the corrected evt_trial
    clear evt_trial
    load('evt.mat');
    
    figure, hold on
    subplot(2,1,1), hold on,
    plot(time(((maxtime*(1/4))*sr:((maxtime*(1/4))+50)*sr)), evt_trial(((maxtime*(1/4))*sr:((maxtime*(1/4))+50)*sr)),'b')
    plot(time(((maxtime*(1/4))*sr:((maxtime*(1/4))+50)*sr)), evt_puff(((maxtime*(1/4))*sr:((maxtime*(1/4))+50)*sr)).*0.5,'r')
    plot(time(((maxtime*(1/4))*sr:((maxtime*(1/4))+50)*sr)), evt_delay(((maxtime*(1/4))*sr:((maxtime*(1/4))+50)*sr)).*0.3,'g')
    legend('trigtrial','puffL','delay')
    xlabel('1min time (s)')
    subplot(2,1,2), hold on,
    plot(time(((maxtime*(3/4))*sr:((maxtime*(3/4))+50)*sr)), evt_trial(((maxtime*(3/4))*sr:((maxtime*(3/4))+50)*sr)),'b')
    plot(time(((maxtime*(3/4))*sr:((maxtime*(3/4))+50)*sr)), evt_puff(((maxtime*(3/4))*sr:((maxtime*(3/4))+50)*sr)).*0.5,'r')
    plot(time(((maxtime*(3/4))*sr:((maxtime*(3/4))+50)*sr)), evt_delay(((maxtime*(3/4))*sr:((maxtime*(3/4))+50)*sr)).*0.3,'g')
    title('after correction')
    saveas(gcf,'evtfig_corrected','png')
    %%
    %%     Manual_time correction
%     evt_puff = evt_puff_L + evt_puff_R;
%     time_2obs = 1300
%     figure, hold on
%     plot(time(time_2obs*sr:(time_2obs+50)*sr), evt_trial(time_2obs*sr:(time_2obs+50)*sr),'b')
%     plot(time(time_2obs*sr:(time_2obs+50)*sr), evt_puff(time_2obs*sr:(time_2obs+50)*sr).*0.5,'r')
%     plot(time(time_2obs*sr:(time_2obs+50)*sr), evt_delay(time_2obs*sr:(time_2obs+50)*sr).*0.3,'g')
%     legend('trigtrial','puffL','delay')
%     xlabel('1min time (s)')
%     title('after correction')
    %         saveas(gcf,'evtfig_correction_MAnual_time_obs','png')
    
end
disp('DONE and DONE!')
