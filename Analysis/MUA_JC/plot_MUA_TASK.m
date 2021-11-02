%% threshold Spikes for MUA

clear all, 

if isempty(ls('S*.mat'))
    dat2mat_JC_Script
end

ChanID = 'S2Ch6'

%% load vectors microV and time
load([ChanID '.mat']); % microV vec
load time.mat; % time vec 
load evt.mat
load info.mat

%% Filter data for spikes
[B,A] = butter(2,250/10000,'high'); % high pass filter (250Hz)
microV_filt = filtfilt(B,A, microV);

%% threshold
% thr = -abs(1.5*std(microV_filt(find(evt_delay)))); % threshold is -3std from filtered signal.
thr= -55 %microV
thr_bool = microV_filt < thr; % valeur de voltage inferieur au threshold
thr_line = thr*ones(1,max(size(time))); % trh line for figure

%% Define spk = 31 points 
% find the total Nb of spk 
Nspk_tot = sum((diff(find(thr_bool))>1)); 

% find all the microV idx that passes the threshold 
spk_thr_idx = find(thr_bool); 

% find the microV idx start and end that passes the threshold 
spk_thr_end_bool = diff(find(thr_bool))>1; %00100000100001  = 3spks, the first one crosses thr by 3 points (001)   
spk_thr_start_bool =  logical([1; spk_thr_end_bool(1:end-1)]);  

spk_thr_end_idx = spk_thr_idx(spk_thr_end_bool);
spk_thr_start_idx = spk_thr_idx(spk_thr_start_bool);

% find the pic for each spk 
for Nspk = 1:Nspk_tot
%     disp(['loop running:  '  num2str(Nspk) ' / ' num2str(Nspk_tot)])
    if spk_thr_start_idx(Nspk)== spk_thr_end_idx(Nspk)
        spk_pic_list_idx(Nspk) = spk_thr_end_idx(Nspk); 
    else
        spk_segment= (microV(spk_thr_start_idx(Nspk):spk_thr_end_idx(Nspk)));
        pic_value = min(spk_segment);
        pic_after_start_idx = min(find(pic_value==spk_segment));
        spk_pic_list_idx(Nspk) = spk_thr_start_idx(Nspk) -1 + pic_after_start_idx;

    end
end
% define a time vector list for all spk based on their pic value
spk_time = time(spk_pic_list_idx); 
spk_pic_value = microV_filt(spk_pic_list_idx);
%%
% create a Matrice with all NSpikes each represented as a 31 point spk vector  
for Nspk = 1:10; % Nspk_tot 
spk_vec_mat(:,Nspk) = microV_filt(spk_pic_list_idx(Nspk)-15:1:spk_pic_list_idx(Nspk)+15);
spk_time_mat(:,Nspk) = time(spk_pic_list_idx(Nspk)-15:1:spk_pic_list_idx(Nspk)+15); 
end
%%


% plot 
%  close all, 
figure, hold on, 
    plot(time,thr_bool*-400,'c','LineWidth',1)
    plot(time, microV_filt,'color',[0.2 0.2 0.2])
    %plot(time,microV,'color',[0.8 0.8 0.8])
    plot(time, thr_line, 'm','LineWidth',1.5)
    plot(spk_time, spk_pic_value, '*k')
    ylim([-200 100])
    xlim([0 50])
% 


%% PLOT THRESH

% reduce plot to a certain time interval
minute_start = 0.0001 
minute_end = 1
SFreq=info.info_freq_parameters.amplifier_sample_rate
sec=SFreq; 
maxtime = max(time)/60; % in minute
disp(['max time is ' num2str(maxtime) 'minutes'])

epoch = [minute_start*sec*60:minute_end*sec*60]; 

% plot
% close all, 
    figure, hold on,
    

    
    plot(time(epoch),thr_bool(epoch)*-400,'c','LineWidth',1)
    plot(time(epoch), microV_filt(epoch),'color',[0.2 0.2 0.2])
    %plot(time,microV,'color',[0.8 0.8 0.8])
    plot(time(epoch), thr_line(epoch), 'm','LineWidth',1.5)
    
    plot(time(epoch), (evt_lick_L(epoch)*-60),'r','LineWidth',1)
    plot(time(epoch), (evt_lick_R(epoch)*-60),'b','LineWidth',1)
    plot(time(epoch), (evt_delay(epoch)*-60), 'g','LineWidth',2)
    plot(time(epoch), (evt_puff_L(epoch)*-30),'r','LineWidth',2)
    plot(time(epoch), (evt_puff_R(epoch)*-30),'b','LineWidth',2)
    plot(spk_time, spk_pic_value, '*k')

    
    plot(time(epoch), zeros(size(time(epoch))), 'k','LineWidth',2)
 
    
 ylabel('microV')
 xlabel('time (s)')
 legend('spk', ['ChanID' ChanID],'treshold','Lick Left','Lick Right', 'Delay', 'Puff L', 'Puff R')
 ylim([-500 100])   
%  xlim([100 300])
    
DDDD
% 
% %% Separate data with Opto Light ON and OFF
% 
% microV_off = microV_filt(digIN==0);
% microV_on = microV_filt(digIN==1);
% 
% %% spk count total
% 
% Nspk_tot = sum((diff(find(thr_bool))>1))
% 
% %% spk count total Ligth vs No light
% thr_bool_lightON = microV_on < thr;
% thr_bool_lightOFF = microV_off < thr;
% 
% NSpk_tot_lightON  = sum((diff(find(thr_bool_lightON))>1))
% NSpk_tot_lightOFF = sum((diff(find(thr_bool_lightOFF))>1))
% 
% % create 1sec bin
% bintime=0.2; % sec
% SampFreq= 20000;
% binsize= SampFreq*bintime; % 1sec bin
% idxON=find(digIN);
% nbin_tot=round(size(idxON,1)/binsize)-2
% countSPk_vec=[];
% 
% for nbin=1:nbin_tot
%     
%     idxbin = 100 + nbin*binsize: nbin*binsize+binsize -100;
%     NSpk_bin_lightON  = sum((diff(find(thr_bool_lightON(idxbin)))>1));
%     NSpk_bin_lightOFF  = sum((diff(find(thr_bool_lightOFF(idxbin)))>1));
%     
%     countSPk_vec(nbin,:) = [NSpk_bin_lightON'; NSpk_bin_lightOFF'];
% end
% 
% %% STATS
% ON  = countSPk_vec(:,1);
% OFF = countSPk_vec(:,2);
% 
% % Wilcoxon/Mann-Withney test (NON PARAMETRIC):
% %if H=1 ==>  "medians are not equal" (two-tailed test, default)
% [P0, H0] = ranksum(ON,OFF)
% 
% %if H=1 ==>  "median of X is greater than median of Y" (right-tailed test)
% [P1, H1] = ranksum(ON,OFF, 'tail', 'right')
% 
% %if H=1 ==>  "median of X is greater than median of Y" (right-tailed test)
% [P2, H2] = ranksum(ON,OFF, 'tail', 'left')
% %%
% 
% ChanSpkCount = [NSpk_tot_lightON, NSpk_tot_lightOFF ];
% 
% SigBarY = [max(ChanSpkCount)+(max(ChanSpkCount)/10) max(ChanSpkCount)+(max(ChanSpkCount)/10)];
% 
% figure, hold on, 
% c = categorical({'LightON','LightOFF'});
% bar(c, ChanSpkCount,0.5,'k')
% plot(SigBarY, 'r')
% ylabel('#spikes')
% ylim(gca,[0 SigBarY(1)+ (SigBarY(1)/5)])
% legend(['Wilcoxon 2tailed H=' num2str(H0)],['P=' num2str(P0)])
% 




