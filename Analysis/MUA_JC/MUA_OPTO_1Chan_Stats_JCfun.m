function MUA_OPTO_1Chan_Stats_JCfun(ChanID, post_pre_task)
% function MUA_OPTO_1Chan_Stats_JCfun(ChanID)
% Find All BINs of Opto ON and compare to the equivalent number of BINs of opto OFF.
% Nbin and Binsize should be the same
% Ttest firing rate ON vs OFF
% Written by JC 11/13/2018
% Last modified by JC 11/14/2018


% ChanID='S1Ch6'

%% loading vectors of interest

load('StimSize_Epoch.mat')
load('time.mat');
load('info.mat'); MouseID = info.info_notes.MouseID; Day = info.info_notes.Day;
load('evt.mat', 'evt_opto', 'evt_trial', 'evt_lick_L', 'evt_lick_R');
evt_lick = evt_lick_L + evt_lick_R; clear evt_lick_R evt_lick_L;
load('Epochs_pre_post_task_st_end.mat', 'idx_preStim_st_end', 'idx_task_st_end','idx_postStim_st_end');

if post_pre_task==1 % PostStim
    evo =   evt_opto(idx_postStim_st_end(1):idx_postStim_st_end(2));
elseif post_pre_task==2 % preStim
    evo =   evt_opto(idx_preStim_st_end(1):idx_preStim_st_end(2));
elseif post_pre_task==3% TaskStim
    evo =   evt_opto(idx_task_st_end(1):idx_task_st_end(2));
end


for datTyp=1:2
    if datTyp==1
        sub_or_raw = 'raw'
    elseif datTyp==2
        sub_or_raw = 'sub'
    end
    
    load([ChanID '_' sub_or_raw '.mat'], 'data','sr');
    %% re-named variables of interest
    
    if post_pre_task == 1
        idx_st_end = idx_postStim_st_end;
        EpochID = 'post'
    elseif post_pre_task == 2
        idx_st_end = idx_preStim_st_end;
        EpochID = 'pre'
    elseif post_pre_task == 3
        idx_st_end = idx_task_st_end;
        EpochID = 'task'
    end
    
    dat =   data(idx_st_end(1):idx_st_end(2));
    evo =   evt_opto(idx_st_end(1):idx_st_end(2));
    tim =   time(idx_st_end(1):idx_st_end(2));
    
    %% Filter data for spikes
    [B,A] = butter(2,750/10000,'high'); % high pass filter (750Hz)
    dat_filt = filtfilt(B,A, dat);
    
    %% threshold
    thr = -round(abs(3*std(dat_filt)))% threshold used -4std from filtered signal.
    % thr = -50 % micron Volts
    thr_bool = dat_filt < thr; % valeur de voltage inferieur au threshold
    thr_line = thr*ones(1,max(size(tim))); % trh line for figure
    
    %% Separate data with Opto Light ON and OFF
    dat_off = dat_filt(evo==0);
    dat_on = dat_filt(evo==1);
    
    limsize = min(size(dat_off,2),size(dat_on,2))
    dat_off = dat_off(1:limsize);% dat_off = dat_off(end-limsize:end);
    dat_on = dat_on(1:limsize);
    
    %% spk count total
    Nspk_tot = sum((diff(find(thr_bool))>1))
    
    %% mean size of a stim pulse
    StimSize = median(diff(find(diff(find(evo))>1)))
    StimTime= roundn((StimSize/sr),-1)
    
    %% spk count total Ligth vs No light
    thr_bool_lightON = dat_on < thr;
    thr_bool_lightOFF = dat_off < thr;
    
    NSpk_tot_lightON  = sum((diff(find(thr_bool_lightON))>1))
    NSpk_tot_lightOFF = sum((diff(find(thr_bool_lightOFF))>1))
    
    % create 1sec bin
    bintime=StimTime % sec
    binsize= sr*bintime;
    
    idxON=find(evo);
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
    
    MeanON = mean(ON)
    MeanOFF = mean(OFF)
    
    stdON = std(ON)
    stdOFF = std(OFF)
    
    % symetric 2 tailed tests:
    [Ht, Pt] = ttest(ON,OFF) % T-test (PARAMETRIC) : if H=0 ==> equals Means
    [Pr, Hr] = ranksum(ON,OFF) % Wilcoxon/Mann-Withney (NON PARAMETRIC): if H=0 ==> equal medians
    
    % Asymetric 1 sided test:
    % if H=1 ==>  "median of X is greater than median of Y" (right-tailed test)
    [PrONoff, HrONoff] = ranksum(ON,OFF, 'tail', 'right')
    %if H=1 ==>  "median of X is less than median of Y" (left-tailed test)
    [PrOFFon, HrOFFon] = ranksum(ON,OFF, 'tail', 'left')
    
    %% Plot
    figure(1), close(1), figure(1),
    
    subplot(1,10,[1:7])
    xlim([tim(1)+45 tim(1)+65])
    ylim([-100 100])
    hold on, plot(tim, dat_filt,'color',[0.2 0.2 0.2])
    hold on,plot(tim, (evo*100)-150, 'c' , 'LineWidth',2)
    hold on,plot(tim, thr_line, '--r', 'LineWidth',1)
    ylabel('microVolts','FontSize', 11);
    xlabel('time (s)')
    legend('VM recorded traces','Laser ON', 'threshold')
    title([MouseID ' ' Day ' ' EpochID ' ' ChanID ' ' sub_or_raw  ], 'FontSize', 12)
    
    subplot(1,10,[9:10])
    barwith = 0.4
    pos=[0.5; 1.5]
    y= [ MeanON/bintime; MeanOFF/bintime]
    std_dev =   [stdON; stdOFF]
    hold on, bar(pos(1), y(1), barwith, 'FaceColor','c','EdgeColor','k','LineWidth',1.5);
    hold on, bar(pos(2), y(2), barwith, 'FaceColor','w','EdgeColor','k','LineWidth',1.5);
    hold on, errorbar(pos(1), y(1), std_dev(1),'.','Color','k','LineWidth',1);
    hold on, errorbar(pos(2), y(2), std_dev(2),'.','Color','k','LineWidth',1);
    set(gca, 'XTick',[0.5 1.5], 'XTickLabel',{'ON','OFF'}); % set(gca,'YLim', [0 50]);
    ylabel('Firing rate (Mean/pusle)','FontSize', 11);
    legend(['Optostim  = pulse ON (' num2str(bintime) 'sec)'] ,['Baseline = pulse OFF(' num2str(bintime) 'sec)'] ,'Location','best')
    title( ['P=' num2str(Pr) ' (rank)'], 'FontSize', 12)
    
    %% SAVING
    
    mkdir('.\MUA')
    saveas(gcf, ['.\MUA\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr' num2str(thr) '_' sub_or_raw '_' EpochID '.png'])
    if HrOFFon ==1
        saveas(gcf, ['fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr' num2str(thr) '_' sub_or_raw '_' EpochID '.png'])
    end
    
    
    if post_pre_task==1
        mkdir('D:\JC_Figures\MUA\PostStim')
        % to save all fig: saveas(gcf, ['D:\JC_Figures\MUA\PostStim\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr) '_' sub_or_raw '_' EpochID '.png'])
        % Save the significant chanel only:
        if HrOFFon ==1
            mkdir('D:\JC_Figures\MUA\PostStim\HrOFFon'); mkdir('D:\JC_Figures\MUA\PostStim\HrOFFon\raw');
            if sub_or_raw == 'sub'
                saveas(gcf, ['D:\JC_Figures\MUA\PostStim\HrOFFon\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr)  '_' sub_or_raw '_' EpochID '.png'])
            elseif sub_or_raw == 'raw'
                saveas(gcf, ['D:\JC_Figures\MUA\PostStim\HrOFFon\raw\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr)  '_' sub_or_raw '_' EpochID '.png'])
            end
        elseif HrONoff ==1
            mkdir('D:\JC_Figures\MUA\PostStim\HrONoff'); mkdir('D:\JC_Figures\MUA\PostStim\HrONoff\raw');
            if sub_or_raw == 'sub'
                saveas(gcf, ['D:\JC_Figures\MUA\PostStim\HrONoff\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr)  '_' sub_or_raw '_' EpochID '.png'])
            elseif sub_or_raw == 'raw'
                saveas(gcf, ['D:\JC_Figures\MUA\PostStim\HrONoff\raw\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr)  '_' sub_or_raw '_' EpochID '.png'])
            end
        end
        
    elseif post_pre_task==2
        mkdir('D:\JC_Figures\MUA\PreStim')
        % to save all fig: saveas(gcf, ['D:\JC_Figures\MUA\PreStim\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr) '_' sub_or_raw '_' EpochID '.png'])
        % Save the significant chanel only:
        if HrOFFon ==1
            mkdir('D:\JC_Figures\MUA\PreStim\HrOFFon'); mkdir('D:\JC_Figures\MUA\PreStim\HrOFFon\raw');
            if sub_or_raw == 'sub'
                saveas(gcf, ['D:\JC_Figures\MUA\PreStim\HrOFFon\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr)  '_' sub_or_raw '_' EpochID '.png'])
            elseif sub_or_raw == 'raw'
                saveas(gcf, ['D:\JC_Figures\MUA\PreStim\HrOFFon\raw\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr)  '_' sub_or_raw '_' EpochID '.png'])
            end
        elseif HrONoff ==1
            mkdir('D:\JC_Figures\MUA\PreStim\HrONoff'); mkdir('D:\JC_Figures\MUA\PreStim\HrONoff\raw');
            if sub_or_raw == 'sub'
                saveas(gcf, ['D:\JC_Figures\MUA\PreStim\HrONoff\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr)  '_' sub_or_raw '_' EpochID '.png'])
            elseif sub_or_raw == 'raw'
                saveas(gcf, ['D:\JC_Figures\MUA\PreStim\HrONoff\raw\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr)  '_' sub_or_raw '_' EpochID '.png'])
            end
        end
        
    elseif post_pre_task==3
        mkdir('D:\JC_Figures\MUA\TaskStim')
        % to save all fig: saveas(gcf, ['D:\JC_Figures\MUA\TaskStim\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr) '_' sub_or_raw '_' EpochID '.png'])
        % Save the significant chanel only:
        if HrOFFon ==1
            mkdir('D:\JC_Figures\MUA\TaskStim\HrOFFon'); mkdir('D:\JC_Figures\MUA\TaskStim\HrOFFon\raw');
            if sub_or_raw == 'sub'
                saveas(gcf, ['D:\JC_Figures\MUA\TaskStim\HrOFFon\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr)  '_' sub_or_raw '_' EpochID '.png'])
            elseif sub_or_raw == 'raw'
                saveas(gcf, ['D:\JC_Figures\MUA\TaskStim\HrOFFon\raw\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr)  '_' sub_or_raw '_' EpochID '.png'])
            end
        elseif HrONoff ==1
            mkdir('D:\JC_Figures\MUA\TaskStim\HrONoff'); mkdir('D:\JC_Figures\MUA\TaskStim\HrONoff\raw');
            if sub_or_raw == 'sub'
                saveas(gcf, ['D:\JC_Figures\MUA\TaskStim\HrONoff\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr)  '_' sub_or_raw '_' EpochID '.png'])
            elseif sub_or_raw == 'raw'
                saveas(gcf, ['D:\JC_Figures\MUA\TaskStim\HrONoff\raw\fig_' MouseID '_' Day '_' ChanID  '_STAT_MUA_OPTO_thr'  num2str(thr)  '_' sub_or_raw '_' EpochID '.png'])
            end
        end
    end
end

