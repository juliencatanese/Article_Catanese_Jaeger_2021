% Plot_Artifact_removal_1ChanRef_JC_Script
% S2Ch7 = ref (the one with the less number of Spikes after AutoClustering)
% Written by Julien Catanese
% Last Updated: 10/08/2018 

clear all
list_nspk= [];
path = pwd;

%% TRY WITH THE MEAN THIS WAY IS NOT GREAT:  
spikes_ChID = dir([path '\*_raw_spikes.mat']);
if ~isempty(spikes_ChID)
    for nsc= 1:max(size(spikes_ChID))
        load (spikes_ChID(nsc).name) ;
        Nspk=size(index,2);
        list_nspk= [list_nspk Nspk];
    end
    
    list_nspk
    min(list_nspk)
    idxref=find(list_nspk==min(list_nspk))
    refID=spikes_ChID(idxref).name
    save('RefChan', 'refID');
end

load ('info.mat')
load('RefChan', 'refID');
time_st= 3; %in minute
time_end=time_st+0.2; % in minute
RefChan=refID(1:5)

raw_ChID = dir([path '\S*Ch*_raw.mat']);
for nrc= 1:max(size(raw_ChID))
    ShChID= raw_ChID(nrc).name
    Sh= ShChID(1:2)
    Ch= ShChID(3:5)
    
    load([RefChan '_raw.mat']);
    ref=data(sr*60*time_st:sr*60*time_end);
    all_ref=data;
    clear data;
    
    load([num2str(Sh) num2str(Ch) '_raw.mat']);
    ch1=data(sr*60*time_st:sr*60*time_end);
    all_ch1=data;
    clear data;
    
    [B,A] = butter(2,600/10000,'high'); % high pass filter
    
    filt_ref = filtfilt(B,A, ref);
    filt_ch1 = filtfilt(B,A, ch1);
    
    load('time.mat');
    t=time(sr*60*time_st:sr*60*time_end);
    clear time;
    
    %% Compare plot
    close all,
    
    figure,
    
    subplot(2,2,1), hold on,
    % plot raw data
    plot(t, ch1)
    plot(t, ref)
    legend(['ch=' ShChID 'raw]'], ['ref= ' RefChan 'raw'], 'Location','best' )
    ylim([-400 300])
    title([ShChID(1:5)])
    
    subplot(2,2,2), hold on,
    % substraction on Raw data.
    ch1_rawsub = ch1 - ref;
    plot(t, ch1)
    plot(t, ch1_rawsub)
    legend('ch1raw', 'ch1raw-refraw', 'Location','best' )
    ylim([-400 300])
    title([info.info_notes.MouseID ' ' info.info_notes.Day ' ' info.info_notes.Depth]) 
    
    subplot(2,2,3), hold on,
    % filter without substraction
    plot(t, filt_ch1)
    % filtering before substraction
    ch1_filtsub = filt_ch1 - filt_ref;
    plot(t, ch1_filtsub)
    % filtering after substraction
    ch1_rawsub_filt = filtfilt(B,A, ch1_rawsub);
    plot(t, ch1_rawsub_filt+200)
    
    ylim([-400 300])
    legend( ' filter without substraction', 'filtering before substraction  ', 'filtering after substraction ' , 'Location','best' )
    
    
    
    % threshold
    STD_THR =2.5;
    all_ch1_rawsub_filt = filtfilt(B,A, (all_ch1-all_ref));
    thr = -abs(STD_THR*std(all_ch1_rawsub_filt)); % threshold is -3std from filtered signal.
    thr_line1 = thr*ones(1,max(size(t))); % trh line for figure
    
    thr = -abs(25*std(ch1_rawsub_filt)); % threshold is -3std from filtered signal.
    thr_line2 = thr*ones(1,max(size(t))); % trh line for figure
    
    % thr_bool = ch1_rawsub_filt < thr; % valeur de voltage inferieur au threshold
    
    subplot(2,2,4), hold on,
    plot(t, ch1_rawsub_filt)
    plot(t, thr_line1)
    plot(t, thr_line2)
    ylim([-400 300])
    title(['thr=' num2str(STD_THR) 'std   to   thr=25std'  ])
    saveas(gcf, ['fig2print_' ShChID(1:5) '_Artifact1_removal_trh'],'jpg')
    
    disp(['SAVING: ' ShChID(1:5) '_sub.mat'])
    
    
    if ShChID(1:5)== RefChan
        data = all_ch1;
        save([ShChID(1:5) '_sub.mat'],'data','sr','unit')
    else
        data = all_ch1-all_ref;
        save([ShChID(1:5) '_sub.mat'],'data','sr','unit')
    end
    disp('saved')
end

%%
%
% sr=20000;
% unit='microVolts'
% raw_ChID = dir([path '\S*Ch*_raw.mat']);
% for nrc= 1:max(size(raw_ChID))
%     ShChID= raw_ChID(nrc).name
%     load([ShChID(1:5) '_sub.mat'])
%     save([ShChID(1:5) '_sub.mat'],'data','sr','unit')
% end
