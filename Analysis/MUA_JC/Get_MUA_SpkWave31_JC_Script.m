% Get_MUA_SpkWave31_JC_Script
% Filter data (High pass 250Hz) 
% Threshold (-55 microV) to detect spk
% Find spk pic 
% Extract 15pts before and 15pts after pic (31pts/spk)
% plot the MUA spks
% Save the spk values and times 
% by Julien Catanese in JaegerLab 
% last update: 06-Apr-2018 

%% threshold Spikes for MUA

clear all, 

if isempty(ls('S*Ch*_raw.mat'))
    dat2mat_JC_Script
end

% for 

ChanID = 'S4Ch1'

%% load vectors microV and time
load([ChanID '_raw.mat']); % microV vec
load time.mat; % time vec 
load evt.mat
load info.mat
microV =data'; 


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
spk1_pic_time = time(spk_pic_list_idx); 
spk1_pic_value = microV_filt(spk_pic_list_idx);
%%
% create a Matrice with all NSpikes each represented as a 31 point spk vector  
spk31_vec=; spk31_time=[]; 
for Nspk = 1: Nspk_tot 
spk31_vec(:,Nspk) = microV_filt(spk_pic_list_idx(Nspk)-15:1:spk_pic_list_idx(Nspk)+15);
spk31_time(:,Nspk) = time(spk_pic_list_idx(Nspk)-15:1:spk_pic_list_idx(Nspk)+15); 
end

%% plot
figure, plot(spk31_vec(:,1:1000))

%% SAVE
save ([ChanID '_spk_MUA_wave31_JC.mat'], 'spk31_vec', 'spk31_time','spk1_pic_value','spk1_pic_time')


 