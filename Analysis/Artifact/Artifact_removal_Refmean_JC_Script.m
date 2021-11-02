% Artifact_removal_Refmean_JC_Script
% S2Ch7 = ref (the one with the less number of Spikes after AutoClustering)

clearvars -except MouseID FolderID nf 
path = pwd; 

%% TRY WITH THE MEAN
% if isempty(dir('*Ref.mat'))

load('time.mat');
load ('info.mat');

if isempty(dir('Ref_mean.mat'))
    MAT=zeros(size(time'));
    rawChID = dir([path '\S*Ch*_raw.mat']);
    for nc= 1:max(size(rawChID))
        disp(['Computing mean REf: ' round(num2str( round((nc/32)*100))) '% complete'])
        clear data;
        load (rawChID(nc).name) ;
        MAT=MAT+data;
    end
    %
    data_ref = MAT/max(size(rawChID));
    unit_ref= 'microV (mean all chan)';
    save('Ref_mean', 'data_ref', 'sr','unit_ref');
else
       disp('Ref_mean.mat already saved') 
end

%%
load('Ref_mean.mat');
all_ref=data_ref;
clear data_ref;

%% loop for each raw Channel 
raw_ChID = dir([path '\S*Ch*_raw.mat']);
for nrc= 1:max(size(raw_ChID))
    ShChID= raw_ChID(nrc).name;
    Sh= ShChID(1:2);
    Ch= ShChID(3:5);
    
    % load raw 
    load([num2str(Sh) num2str(Ch) '_raw.mat'],'data','sr','unit');
    all_ch1 = data;
    clear data;        
    
    % compute and save sub 
    data = all_ch1 - all_ref;    
    save([ShChID(1:5) '_sub.mat'],'data','sr','unit')
    clear data;
    disp(' .sub saved')
   
    %% plotting  
    %  Define beginning and END for plotting a sample  
    time_st= 5; %in minute
    time_end=time_st+0.1; % in minute
    t=time(sr*60*time_st:sr*60*time_end);
    ch1=all_ch1(sr*60*time_st:sr*60*time_end);
    ref=all_ref(sr*60*time_st:sr*60*time_end);
    
    % high pass filter  
    [B,A] = butter(2,600/10000,'high'); 
    
    filt_ref = filtfilt(B,A, ref);
    filt_ch1 = filtfilt(B,A, ch1);
    
    %% Compare plot
    close all,
    figure,
    subplot(2,2,1), hold on,
    % plot raw data
    plot(t, ch1)
    plot(t, ref)
    legend(['ch=' ShChID 'raw]'], ['ref= Mean all raw'], 'Location','best' )
    ylim([-400 200])
    title([ShChID(1:5)])
    
    subplot(2,2,2), hold on,
    % substraction on Raw data.
    ch1_rawsub = ch1 - ref;
    plot(t, ch1)
    plot(t, ch1_rawsub)
    legend('ch1raw', 'ch1raw-refraw', 'Location','best' )
    ylim([-400 200])
    title([info.info_notes.MouseID ' ' info.info_notes.Day ' ' info.info_notes.Depth])
    
    subplot(2,2,3), hold on,
    % filter without substraction
    plot(t, filt_ch1)
    % filtering before substraction
    ch1_filtsub = filt_ch1 - filt_ref;
    plot(t, ch1_filtsub)
    % filtering after substraction
    ch1_rawsub_filt = filtfilt(B,A, ch1_rawsub);
    plot(t, ch1_rawsub_filt)
    
    ylim([-400 200])
    legend( ' filter without substraction', 'filtering before substraction  ', 'filtering after substraction ' , 'Location','best' )
    
    % threshold
    STD_THR =4;
    all_ch1_rawsub_filt = filtfilt(B,A, (all_ch1-all_ref));
    thr = -abs(STD_THR*std(all_ch1_rawsub_filt)); % threshold is -3std from filtered signal.
    thr_line1 = thr*ones(1,max(size(t))); % trh line for figure
    
    thr = -abs(16*std(ch1_rawsub_filt)); % threshold is -3std from filtered signal.
    thr_line2 = thr*ones(1,max(size(t))); % trh line for figure
    
    % thr_bool = ch1_rawsub_filt < thr; % valeur de voltage inferieur au threshold
    
    subplot(2,2,4), hold on,
    plot(t, ch1_rawsub_filt)
    plot(t, thr_line1)
    plot(t, thr_line2)
    ylim([-400 200])
    title(['thr=' num2str(STD_THR) 'std   to   thr=20std'  ])
    legend('filtering after substraction')
    
    saveas(gcf, ['fig2print_' ShChID(1:5) '_Artifact1_removal_v2'],'jpg')
    
    disp(['SAVING: ' ShChID(1:5) '_sub.mat'])    
    disp('figure saved')
end

