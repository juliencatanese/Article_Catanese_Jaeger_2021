%% threshold Spikes for MUA

clear all, close all;

if isempty(ls('*.mat'))
    Get_data_per_chan_JC_Script
end

clear all, close all;

%% load vectors microV and time
load S1Ch6_sub.mat; % microV vec
load time.mat; % time vec 
load evt.mat
microV=data; 
digIN=evt_opto;

%% Optional: Manually reduce size of microV and time
microV=data; 
timeori=time; 
microV=microV(find(time==400):find(time==600));
digIN=digIN(find(time==400):find(time==600));
time=time(find(time==400):find(time==600));

%% Filter data for spikes
[B,A] = butter(2,750/10000,'high'); % high pass filter (750Hz)
microV_filt = filtfilt(B,A, microV);

%% threshold
thr = -abs(3*std(microV_filt)); % threshold is -3std from filtered signal.
thr_bool = microV_filt < thr; % valeur de voltage inferieur au threshold
thr_line = thr*ones(1,max(size(time))); % trh line for figure

%% PLOT THRESH
close all, 
    figure, hold on,
%     plot(time,microV,'color',[0.8 0.8 0.8])
%     plot(time,thr_bool*-400,'c')
    plot(time, microV_filt,'color',[0.2 0.2 0.2])
%     plot(time, thr_line, 'r')
    plot(time, (digIN*100)-50)

 ylabel('microV')
 xlabel('time (s)')
%  legend('SNr electrode','OPTO SNr ARCH (ON>0 vs OFF<0)')
legend('VM electrode','OPTO ChR2 (ON>0 vs OFF<0)')

    
%% Separate data with Opto Light ON and OFF

microV_off = microV_filt(digIN==0);
microV_on = microV_filt(digIN==1);

%% spk count total

Nspk_tot = sum((diff(find(thr_bool))>1))

%% spk count total Ligth vs No light
thr_bool_lightON = microV_on < thr;
thr_bool_lightOFF = microV_off < thr;

NSpk_tot_lightON  = sum((diff(find(thr_bool_lightON))>1))
NSpk_tot_lightOFF = sum((diff(find(thr_bool_lightOFF))>1))

% create 1sec bin
bintime=0.2; % sec
SampFreq= 20000;
binsize= SampFreq*bintime; % 1sec bin
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

figure, hold on, 
c = categorical({'LightON','LightOFF'});
bar(c, ChanSpkCount,0.5,'k')
plot(SigBarY, 'r')
ylabel('#spikes')
ylim(gca,[0 SigBarY(1)+ (SigBarY(1)/5)])
legend(['Wilcoxon 2tailed H=' num2str(H0)],['P=' num2str(P0)])





