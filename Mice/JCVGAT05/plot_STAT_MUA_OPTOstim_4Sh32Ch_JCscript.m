%% threshold Spikes for MUA

clear all, close all;

%%
plotRawTraces = 0
plotbins = 0
plotspktime = 0;
plotRASTER=0;

%% Define Data Folders and Figure Folder.  
DataFolder = 'D:\JC_Data\Acute\Awake\JC-VGAT05\day2\'
FigureFolder = 'C:\Users\JCATANE\Documents\MATLAB'

%% create List Folder to process 
P = dir(DataFolder)
List_Folder={};
for nN = 1:max(size(P))
    if P(nN).isdir == 1
        List_Folder = [List_Folder P(nN).name];
    end
end
List_Folder = List_Folder(3:end)

%% LOOP ACROSS FOLDER
for nF = 1:max(size(List_Folder))
    Folder= [ P(1).folder '\' List_Folder{nF}]
    
    %% load channel data matrices (chan_data_mat)
    load ([Folder '\AllChan.mat'], 'Chan_data_mat', 'Chan_ID')
    
    %% Create vectors
    time = Chan_data_mat(:,end);
    digin2 = Chan_data_mat(:,end-1);
    for ch=1:32; pause(0.01);
        chan1 = Chan_data_mat(:,ch);
        Chan_ID(:,ch);
        
        %% Filter data for spikes
        [B,A] = butter(2,250/10000,'high'); % high pass filter
        chan1_filt = filtfilt(B,A,chan1);
        
        %% threshold
        thr = -abs(3*std(chan1_filt)); % threshold is -3std from filtered signal.
        thr_bool = chan1_filt < thr; % valeur de voltage inferieur au threshold
        thr_line = thr*ones(1,max(size(time))); % trh line for figure
        
        %% PLOT THRESH
        if plotRawTraces
            figure, hold on,
            plot(time,chan1,'color',[0.8 0.8 0.8])
            plot(time,thr_bool*-400,'c')
            plot(time, chan1_filt,'color',[0.2 0.2 0.2])
            plot(time, thr_line, 'r')
            plot(time, digin2*100) % opto stim
        end
        %% Separate data with Opto Light ON and OFF
        chan1_off = chan1_filt(digin2==0);
        chan1_on = chan1_filt(digin2==1);
       
        
        %% spk count total Ligth vs No light
        thr_bool_lightON = chan1_on < thr;
        thr_bool_lightOFF = chan1_off < thr;
        
                %% spk count total
        
        Nspk_tot = sum((diff(find(thr_bool))>1))
        NSpk_tot_lightOFF  = sum(diff(thr_bool_lightOFF)<0)
        NSpk_tot_lightON  = sum(diff(thr_bool_lightON)<0)
       
        
        %% Define bins length (in sec)
        bintime=0.2; % (sec)
        SampFreq= 20000; %(Hz) % obtain by Get_Data_Allrhd_JC(Folder)
        binsize= SampFreq*bintime; % Xsec bin
        idxON=find(digin2);
        nbin_tot=round(size(idxON,1)/binsize)-2
        countSPk_vec=[];
        
        
        %% find idx of spk during Light ON (idx in vector A).
        A=find(digin2==1); %= chan1_filt(digin2==1); = opto stim ON ___---___---___---___
        B=find(thr_bool_lightON) ;%thr_bool_lightON= chan1_on < thr; = All points that pass thr when stim ON -!-!!!---!!--!!!--
        C=A(B);
        C1=zeros(size(time)); C2=C1;
        C1(C)=1;
        
        %% Spk become a single idx (middle of the segment)
        binspkstart=[]; binspkend=[];
        binspkstart= find(diff(thr_bool_lightON)>0); % diff(X), for a vector X, is [X(2)-X(1)  X(3)-X(2) ... X(n)-X(n-1)]
        binspkend=find(diff(thr_bool_lightON)<0);
        % in case the size mismatch beetween binspkstart and binspkend; 
        if thr_bool_lightON(end) ==1  
            binspkstart(end)=[]; 
        end
        if  thr_bool_lightON(1) ==1
            binspkend(1)=[];   %Test= diff([binspkstart(1:5) binspkend(1:5)],1,2); 
        end
        MM= round(mean([binspkstart  binspkend],2)); size(MM)
        MMM=A(MM); %MMM = idx spk (=1 time) in A -!-!-!---!-!--!-!-!--
        C2(MMM)=1;
        
        % plot spk times
        if plotspktime;
            figure,  hold on,
            plot(time, -C1,'r');
            plot(time, -C2*0.5,'k')
        end
        
        %% PLOT RASTER and create a var to be SAVED for later STATs  
        D= []; spk = [];
        if plotRASTER;
            X=1:binsize;
            figure,
            nn=0;
        end
        % loop for each trial (Bin) 
        for nbin=1:nbin_tot; pause(0.01);
            
            % generate bins idx 1:20000 then 20001:40000 ...
            idxbin = 1 + (nbin*binsize): (nbin*binsize)+binsize ;
            
            % count spk in each bin
            NSpk_bin_lightON  = sum((diff(find(thr_bool_lightON(idxbin)))>1));
            NSpk_bin_lightOFF  = sum((diff(find(thr_bool_lightOFF(idxbin)))>1));
            
            % STAT var to be SAVED = spk count/bin (2 vectors ON and OFF);
            countSPk_vec(nbin,:) = [NSpk_bin_lightON'; NSpk_bin_lightOFF'];
            
            % create D for later visualize binning
            if plotbins
                disp ([ num2str(nbin) '/' num2str(nbin_tot) ]);
                D=[D A(idxbin)];
            end
            
            % PLOT RASTER 
            if plotRASTER;
                nn=nn+1;
                [idx_MMM idx_idxbin]=find(MMM==idxbin);
                %                 Ntrial = ones(size(idx_idxbin))* nbin;
                %                 spk = [spk; idx_idxbin Ntrial];

                Y=[]; Y= NaN(size(X));
                Y(idx_idxbin)=nbin;        
                if nn>50
                    hold on, plot(X,Y,'b.'),%plot(X,Y,'cd');
                else    
                    hold on, plot(X,Y,'r.'),%plot(X,Y,'cd');
                end
            end
        end
        
        
  
        
        
        %%  plot for testing binning
        if plotbins
            figure, hold on,
            plot(time, digin2*100,'c');
            plot(time, -C1*250,'r');
            plot(time, chan1_filt,'color',[0.2 0.2 0.2]);
            C2(D)=1;
            plot(time, C2*150,'g');
        end
        
        %         figure, hold on
        %         for nbin=1:nbin_tot
        %             plot (idxON(100 + nbin*binsize: nbin*binsize+binsize -100)); % change bin time to see the changes in plot
        %         end
        
        %% STATS
        ON  = countSPk_vec(:,1);
        OFF = countSPk_vec(:,2);
        
        % Wilcoxon/Mann-Withney test (NON PARAMETRIC):
        %if H=1 ==>  "medians are not equal" (two-tailed test, default)
        [P0, H0] = ranksum(ON,OFF);
        
        %if H=1 ==>  "median of X is greater than median of Y" (right-tailed test)
        [P1, H1] = ranksum(ON,OFF, 'tail', 'right');
        
        %if H=1 ==>  "median of X is greater than median of Y" (right-tailed test)
        [P2, H2] = ranksum(ON,OFF, 'tail', 'left');
        
        disp(['(ON>OFF) H1=' num2str(H1) '   vs    (OFF>ON) H2=' num2str(H2)])
        
        % TTEST (PARAMETRIC)
        %if H=1 ==> the null hypothesis (equal means) can be rejected at the 5% level.
        % [H,P] = ttest(ON, OFF)
        % % if H=1 ==>
        % [Hks,Pks] = kstest(ON-OFF/mean(ON-OFF))
        % figure, histogram((ON-OFF)/(mean(ON-OFF)))
        %
        % [Hks,Pks] = kstest(ON/mean(ON))
        % figure, histogram((ON)/(mean(ON)))
        %
        % [Hks,Pks] = kstest(OFF/mean(OFF))
        % figure, histogram((OFF)/(mean(OFF)))
        
        ChanSTAT_H_P(ch,:) = [H0 P0];
        
        Chan_mean_bin_SpkCount(ch,:) = mean(countSPk_vec,1);
        
        Chan_std_bin_SpkCount(ch,:) = std(countSPk_vec,1);
        
        Chan_tot_bin_SpkCount(ch,:) = sum(countSPk_vec,1);
        
        Chan_tot_SpkCount(ch,:) = [Nspk_tot, NSpk_tot_lightON, NSpk_tot_lightOFF ];
        
        Chan_Norm = Chan_mean_bin_SpkCount(:,:)./Chan_mean_bin_SpkCount(:,2);
        
        
        
    end
    
    %% SAVE stats var
    save([Folder '/STATS_LightONvsOFF.m'],'ChanSTAT_H_P','countSPk_vec','Chan_mean_bin_SpkCount','Chan_std_bin_SpkCount','Chan_tot_SpkCount','Chan_tot_bin_SpkCount')
    
    %% PLOT stats and SAVE fig
    
    %% Plot relative to OFF Mean SPk Count per Bin + STAT
    colormap(cool),
    Fig1= figure,
    SigBarX = [1 1]*max(max(Chan_Norm))*1.1;
    
    
    subplot(1,4,1), B1=barh(-[1 2 3 4 5 6 7 8],Chan_Norm(1:8,:) ,0.9, 'grouped')
    
    SigBarY1 = -[find(ChanSTAT_H_P(1:8,1))-0.15 find(ChanSTAT_H_P(1:8,1))+0.15]
    hold on, plot(SigBarX,SigBarY1, '-k', 'LineWidth',2)
    hold on, plot(SigBarX(1)*1.05, mean(SigBarY1,2), '*k')
    
    subplot(1,4,2), barh(-[1 2 3 4 5 6 7 8],Chan_Norm(9:16,:) ,0.9, 'grouped')
    
    SigBarY2 = -[find(ChanSTAT_H_P(9:16,1))-0.15 find(ChanSTAT_H_P(9:16,1))+0.15]
    hold on, plot(SigBarX,SigBarY2, '-k', 'LineWidth',2)
    hold on, plot(SigBarX(1)*1.05, mean(SigBarY2,2), '*k')
    
    subplot(1,4,3), barh(-[1 2 3 4 5 6 7 8],Chan_Norm(17:24,:),0.9, 'grouped')
    
    SigBarY3 = -[find(ChanSTAT_H_P(17:24,1))-0.15 find(ChanSTAT_H_P(17:24,1))+0.15]
    hold on, plot(SigBarX,SigBarY3, '-k', 'LineWidth',2)
    hold on, plot(SigBarX(1)*1.05, mean(SigBarY3,2), '*k')
    
    subplot(1,4,4), barh(-[1 2 3 4 5 6 7 8],Chan_Norm(25:32,:),0.9, 'grouped')
    
    SigBarY4 = -[find(ChanSTAT_H_P(25:32,1))-0.15 find(ChanSTAT_H_P(25:32,1))+0.15]
    hold on, plot(SigBarX,SigBarY4, '-k', 'LineWidth',2)
    hold on, plot(SigBarX(1)*1.05, mean(SigBarY4,2), '*k')
     
    subplot(1,4,2),
    title([Folder(27:32) ' ' Folder(34:37) ' STAT OPTO ' List_Folder{nF} ' bin' num2str(bintime) 's Sh4Ch32'])
    
    hold off,
    legend('ON', 'OF');
    legend = legend(gca);
    set(legend,'Position',[0.858271329983743 0.902978388501941 0.123214284756354 0.0869047596341086]);
    clear legend;
    
    %%
    saveas(Fig1, [Folder '/STAT_LightONvsOFF_Probe.jpeg'])
    saveas(Fig1, [ FigureFolder '/Figures/' Folder(27:32) '_' Folder(34:37) '_STAT_OPTO_' List_Folder{nF} '_bin' num2str(bintime) 's_Sh4Ch32_grouped.fig'])
    saveas(Fig1, [ FigureFolder '/Figures/' Folder(27:32) '_' Folder(34:37) '_STAT_OPTO_' List_Folder{nF} '_bin' num2str(bintime) 's_Sh4Ch32_grouped.jpeg'])
    
    
    
    %% Plot Mean SPk Count per Bin + STAT
    colormap(cool),
    Fig2= figure,
    SigBarX = [1 1]*max(max(Chan_mean_bin_SpkCount))*1.1;
    
    
    subplot(1,4,1), B1=barh(-[1 2 3 4 5 6 7 8],Chan_mean_bin_SpkCount(1:8,:) ,0.9, 'grouped')
    
    SigBarY1 = -[find(ChanSTAT_H_P(1:8,1))-0.15 find(ChanSTAT_H_P(1:8,1))+0.15]
    hold on, plot(SigBarX,SigBarY1, '-k', 'LineWidth',2)
    hold on, plot(SigBarX(1)*1.05, mean(SigBarY1,2), '*k')
    
    subplot(1,4,2), barh(-[1 2 3 4 5 6 7 8],Chan_mean_bin_SpkCount(9:16,:) ,0.9, 'grouped')
    
    SigBarY2 = -[find(ChanSTAT_H_P(9:16,1))-0.15 find(ChanSTAT_H_P(9:16,1))+0.15]
    hold on, plot(SigBarX,SigBarY2, '-k', 'LineWidth',2)
    hold on, plot(SigBarX(1)*1.05, mean(SigBarY2,2), '*k')
    
    subplot(1,4,3), barh(-[1 2 3 4 5 6 7 8],Chan_mean_bin_SpkCount(17:24,:),0.9, 'grouped')
    
    SigBarY3 = -[find(ChanSTAT_H_P(17:24,1))-0.15 find(ChanSTAT_H_P(17:24,1))+0.15]
    hold on, plot(SigBarX,SigBarY3, '-k', 'LineWidth',2)
    hold on, plot(SigBarX(1)*1.05, mean(SigBarY3,2), '*k')
    
    subplot(1,4,4), barh(-[1 2 3 4 5 6 7 8],Chan_mean_bin_SpkCount(25:32,:),0.9, 'grouped')
    
    SigBarY4 = -[find(ChanSTAT_H_P(25:32,1))-0.15 find(ChanSTAT_H_P(25:32,1))+0.15]
    hold on, plot(SigBarX,SigBarY4, '-k', 'LineWidth',2)
    hold on, plot(SigBarX(1)*1.05, mean(SigBarY4,2), '*k')
    
    
    subplot(1,4,2),
    title([Folder(27:32) ' ' Folder(34:37) ' STAT OPTO ' List_Folder{nF} ' bin' num2str(bintime) 's Sh4Ch32'])
    
    hold off,
    legend('ON', 'OF');
    legend = legend(gca);
    set(legend,'Position',[0.858271329983743 0.902978388501941 0.123214284756354 0.0869047596341086]);
    clear legend;
    
    %
    saveas(Fig2, [Folder '/STAT_LightONvsOFF_Probe.jpeg'])
    saveas(Fig2, [ FigureFolder '/Figures/' Folder(27:32) '_' Folder(34:37) '_STAT_OPTO_' List_Folder{nF} '_bin' num2str(bintime) 's_Sh4Ch32_grouped.fig'])
    saveas(Fig2, [ FigureFolder '/Figures/' Folder(27:32) '_' Folder(34:37) '_STAT_OPTO_' List_Folder{nF} '_bin' num2str(bintime) 's_Sh4Ch32_grouped.jpeg'])
    
    
     %% PLOT Stacked stats Normalized to OFF 
    Fig3= figure, colormap(cool),
    SigBarX = [1 1]*max(max(Chan_mean_bin_SpkCount))*1.1;
    
    subplot(1,4,1), B1=barh(-[1 2 3 4 5 6 7 8],Chan_Norm(1:8,:) ,0.9, 'stacked')
    
    SigBarY1 = -[find(ChanSTAT_H_P(1:8,1))-0.15 find(ChanSTAT_H_P(1:8,1))+0.15]
    %     hold on, plot(SigBarX,SigBarY1, '-k', 'LineWidth',2)
    hold on, plot(sum(max(Chan_Norm(1:8,:)))*1.15, mean(SigBarY1,2), '*k')
    
    subplot(1,4,2), barh(-[1 2 3 4 5 6 7 8],Chan_Norm(9:16,:) ,0.9, 'stacked')
    
    SigBarY2 = -[find(ChanSTAT_H_P(9:16,1))-0.15 find(ChanSTAT_H_P(9:16,1))+0.15]
    %     hold on, plot(SigBarX,SigBarY2, '-k', 'LineWidth',2)
    hold on, plot(sum(max(Chan_Norm(9:16,:)))*1.15, mean(SigBarY2,2), '*k')
    
    subplot(1,4,3), barh(-[1 2 3 4 5 6 7 8],Chan_Norm(17:24,:),0.9, 'stacked')
    
    SigBarY3 = -[find(ChanSTAT_H_P(17:24,1))-0.15 find(ChanSTAT_H_P(17:24,1))+0.15]
    %     hold on, plot(SigBarX,SigBarY3, '-k', 'LineWidth',2)
    hold on, plot(sum(max(Chan_Norm(17:24,:)))*1.15, mean(SigBarY3,2), '*k')
    
    subplot(1,4,4), barh(-[1 2 3 4 5 6 7 8],Chan_Norm(25:32,:),0.9, 'stacked')
    
    SigBarY4 = -[find(ChanSTAT_H_P(25:32,1))-0.15 find(ChanSTAT_H_P(25:32,1))+0.15]
    %     hold on, plot(SigBarX,SigBarY4, '-k', 'LineWidth',2)
    hold on, plot(sum(max(Chan_Norm(25:32,:)))*1.15, mean(SigBarY4,2), '*k')
    
    subplot(1,4,2), title([Folder(27:32) ' ' Folder(34:37) ' STAT OPTO ' List_Folder{nF} ' bin' num2str(bintime) 's Sh4Ch32'])
    
    %
    saveas(Fig3, [Folder '/STAT_LightONvsOFF_Probe.jpeg'])
    saveas(Fig3, [ FigureFolder '/Figures/' Folder(27:32) '_' Folder(34:37) '_STAT_OPTO_' List_Folder{nF} '_mean_bin' num2str(bintime) 's_Sh4Ch32_stacked_Normalize.fig'])
    saveas(Fig3, [ FigureFolder '/Figures/' Folder(27:32) '_' Folder(34:37) '_STAT_OPTO_' List_Folder{nF} '_mean_bin' num2str(bintime) 's_Sh4Ch32_stacked_Normalize.jpeg'])
    
    
    %% PLOT Stacked stats Mean Spk Count per bin 
%     Fig4= figure, colormap(cool),
%     SigBarX = [1 1]*max(max(Chan_mean_bin_SpkCount))*1.1;
%     
%     subplot(1,4,1), B1=barh(-[1 2 3 4 5 6 7 8],Chan_mean_bin_SpkCount(1:8,:) ,0.9, 'stacked')
%     
%     SigBarY1 = -[find(ChanSTAT_H_P(1:8,1))-0.15 find(ChanSTAT_H_P(1:8,1))+0.15]
%     %     hold on, plot(SigBarX,SigBarY1, '-k', 'LineWidth',2)
%     hold on, plot(sum(max(Chan_mean_bin_SpkCount(1:8,:)))*1.15, mean(SigBarY1,2), '*k')
%     
%     subplot(1,4,2), barh(-[1 2 3 4 5 6 7 8],Chan_mean_bin_SpkCount(9:16,:) ,0.9, 'stacked')
%     
%     SigBarY2 = -[find(ChanSTAT_H_P(9:16,1))-0.15 find(ChanSTAT_H_P(9:16,1))+0.15]
%     %     hold on, plot(SigBarX,SigBarY2, '-k', 'LineWidth',2)
%     hold on, plot(sum(max(Chan_mean_bin_SpkCount(9:16,:)))*1.15, mean(SigBarY2,2), '*k')
%     
%     subplot(1,4,3), barh(-[1 2 3 4 5 6 7 8],Chan_mean_bin_SpkCount(17:24,:),0.9, 'stacked')
%     
%     SigBarY3 = -[find(ChanSTAT_H_P(17:24,1))-0.15 find(ChanSTAT_H_P(17:24,1))+0.15]
%     %     hold on, plot(SigBarX,SigBarY3, '-k', 'LineWidth',2)
%     hold on, plot(sum(max(Chan_mean_bin_SpkCount(17:24,:)))*1.15, mean(SigBarY3,2), '*k')
%     
%     subplot(1,4,4), barh(-[1 2 3 4 5 6 7 8],Chan_mean_bin_SpkCount(25:32,:),0.9, 'stacked')
%     
%     SigBarY4 = -[find(ChanSTAT_H_P(25:32,1))-0.15 find(ChanSTAT_H_P(25:32,1))+0.15]
%     %     hold on, plot(SigBarX,SigBarY4, '-k', 'LineWidth',2)
%     hold on, plot(sum(max(Chan_mean_bin_SpkCount(25:32,:)))*1.15, mean(SigBarY4,2), '*k')
%     
%     subplot(1,4,2), title([Folder(27:32) ' ' Folder(34:37) ' STAT OPTO ' List_Folder{nF} ' bin' num2str(bintime) 's Sh4Ch32'])
%     saveas(Fig4, [Folder '/STAT_LightONvsOFF_Probe.jpeg'])
%     saveas(Fig4, [ FigureFolder '/Figures/' Folder(27:32) '_' Folder(34:37) '_STAT_OPTO_' List_Folder{nF} '_mean_bin' num2str(bintime) 's_Sh4Ch32_stacked.fig'])
%     saveas(Fig4, [ FigureFolder '/Figures/' Folder(27:32) '_' Folder(34:37) '_STAT_OPTO_' List_Folder{nF} '_mean_bin' num2str(bintime) 's_Sh4Ch32_stacked.jpeg'])
%     
    
    
    
    
    %%
    
    
    
end
List_Folder