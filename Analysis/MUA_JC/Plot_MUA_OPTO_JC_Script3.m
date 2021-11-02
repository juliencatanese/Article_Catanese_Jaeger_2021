% Find All BINs of Opto ON and compare to the equivalent number of BINs of opto OFF.
% Nbin and Binsize should be the same
% Ttest firing rate ON vs OFF
% Written by JC 11/13/2018

%% Define Session list
clear all
close all
cd('D:\JC_Analysis');
SessList = dir(['**/*' '_' 'taskopto' '_' '*']);
NSess= max(size(SessList)) % Number of Sessions
Nol=[]; Nolo=[]; cCR=[]; cCRo=[]; cCL=[]; cCLo=[];

%% loop trhough all Sessions named "taskopto"
for of=1:NSess
    SessID= SessList(of).name
    SessPath=SessList(of).folder;
    cd([SessPath '\' SessID])
    
    if isempty(ls('*.mat'))
        disp('Now doing: Get_data_per_chan_JC_Script')
        Get_data_per_chan_JC_Script ; disp('Done');
    end
    
    
    ChID_list = dir([path '\S*Ch*_raw.mat']);
    for nc= 1:max(size(ChID_list))
        ChanID=ChID_list(nc)
        %% load vectors microV and time
        load([ChanID '_sub.mat'], 'data','sr'); microV=data;
        load('time.mat', 'time');time2=time;
        load('info.mat')
        load('TaskEpoch.mat','taskSTART_idx', 'taskEND_idx', 'optoPostSTART_idx', 'optoPostEND_idx')
        %%
        load('evt.mat');
        digIN=evt_opto;
        evt_lick = evt_lick_L + evt_lick_R;
        evt_rwd=evt_rwd_L + evt_rwd_R;
        evt_puff=evt_puff_L + evt_puff_R;
        
        %% plot Opto Stim in time
        close(1)
        figure(1), plot(time, evt_opto, 'c')
        figure(1), hold on, plot(time, evt_lick, 'r')
        % figure(1), hold on, plot(time, evt_delay, 'g')
        %% define task Start and Task End
        taskEND = task
        
        
        
        %% Optional: Manually reduce size of microV and time
        
        microV=microV(find(time2==time2(end-(800*sr))):find(time2==time2(end-(400*sr))));
        digIN=digIN(find(time2==time2(end-(800*sr))):find(time2==time2(end-(400*sr))));
        time2=time2(find(time2==time2(end-(800*sr))):find(time2==time2(end-(400*sr))));
        
        %% Filter data for spikes
        [B,A] = butter(2,750/10000,'high'); % high pass filter (750Hz)
        microV_filt = filtfilt(B,A, microV);
        
        %% threshold
        thr = -abs(3*std(microV_filt)); % threshold is -3std from filtered signal.
        thr_bool = microV_filt < thr; % valeur de voltage inferieur au threshold
        thr_line = thr*ones(1,max(size(time2))); % trh line for figure
        
        %% PLOT THRESH
        close all,
        figure(2), hold on,
        %     plot(time2,microV,'color',[0.8 0.8 0.8])
        %     plot(time2,thr_bool*-400,'c')
        plot(time2, microV_filt,'color',[0.2 0.2 0.2])
        %     plot(time2, thr_line, 'r')
        plot(time2, (digIN*100)-50)
        
        ylabel('microV')
        xlabel('time2 (s)')
        %  legend('SNr electrode','OPTO SNr ARCH (ON>0 vs OFF<0)')
        legend('VM electrode','OPTO ChR2 (ON>0 vs OFF<0)')
        
        %% Separate data with Opto Light ON and OFF
        
        microV_off = microV_filt(digIN==0);
        microV_on = microV_filt(digIN==1);
        
        limsize = min(size(microV_off,2),size(microV_on,2))
        % microV_off = microV_off(end-limsize:end);
        microV_off = microV_off(1:limsize);
        microV_on = microV_on(1:limsize);
        
        %% spk count total
        
        Nspk_tot = sum((diff(find(thr_bool))>1))
        
        %% mean size of a stim pulse
        StimSize = median(diff(find(diff(find(digIN))>1)))
        StimTime= roundn((StimSize/sr),-1)
        
        %% spk count total Ligth vs No light
        thr_bool_lightON = microV_on < thr;
        thr_bool_lightOFF = microV_off < thr;
        
        NSpk_tot_lightON  = sum((diff(find(thr_bool_lightON))>1))
        NSpk_tot_lightOFF = sum((diff(find(thr_bool_lightOFF))>1))
        
        % create 1sec bin
        bintime=round(StimTime); % sec
        SampFreq= 20000;
        binsize= SampFreq*bintime;
        idxON=find(digIN);
        nbin_tot=round(size(idxON,1)/binsize)-2
        countSPk_vec=[];
        
        for nbin=1:nbin_tot
            
            idxbin = 100 + nbin*binsize: nbin*binsize+binsize -100;
            NSpk_bin_lightON  = sum((diff(find(thr_bool_lightON(idxbin)))>1));
            NSpk_bin_lightOFF  = sum((diff(find(thr_bool_lightOFF(idxbin)))>1));
            
            countSPk_vec(nbin,:) = [NSpk_bin_lightON'; NSpk_bin_lightOFF'];
        end
        
        %% STATS
        ON  = countSPk_vec(:,1);
        OFF = countSPk_vec(:,2);
        
        % Wilcoxon/Mann-Withney test (NON PARAMETRIC):
        %if H=1 ==>  "medians are not equal" (two-tailed test, default)
        [P0, H0] = ranksum(ON,OFF)
        
        %if H=1 ==>  "median of X is greater than median of Y" (right-tailed test)
        [P1, H1] = ranksum(ON,OFF, 'tail', 'right')
        
        %if H=1 ==>  "median of X is greater than median of Y" (right-tailed test)
        [P2, H2] = ranksum(ON,OFF, 'tail', 'left')
        %%
        
        ChanSpkCount = [NSpk_tot_lightON, NSpk_tot_lightOFF ];
        
        SigBarY = [max(ChanSpkCount)+(max(ChanSpkCount)/10) max(ChanSpkCount)+(max(ChanSpkCount)/10)];
        
        figure(3), hold on,
        c = categorical({'LightON','LightOFF'});
        bar(c, ChanSpkCount,0.5,'k')
        plot(SigBarY, 'r')
        ylabel('#spikes', 'FontSize',16)
        ylim(gca,[0 SigBarY(1)+ (SigBarY(1)/5)])
        legend(['Wilcoxon 2tailed H=' num2str(H0)],['P=' num2str(P0)])
        
        
    end
end

        
        
