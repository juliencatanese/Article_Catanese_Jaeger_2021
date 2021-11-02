% Plot_DLCLick_vs_TTLSensor_align
% JC Script 10/02/2020
% cd (DATA_folder)
% load('Ntrial_type.mat')
% load ('evt.mat','evt_trial', 'evt_lick_L', 'evt_lick_R'); load ('time.mat');

%% TO compare TTL Lick Sensor vs DLC Lick detection
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);
TTL_lenght_sec=trig_end(2)/20000 - trig_st(2)/20000
%
idx_trial_start = trig_st(1:end-1);
idx_trial_end = idx_trial_start + 4*20000;
evt_lick = evt_lick_L + evt_lick_R;
%% Get Licks Sensor Evt in 25Hz frames
tr=5
idx_tr = [idx_trial_start(tr):1:idx_trial_end(tr)]
Xt =time(idx_tr); Xt(1)
figure, hold on,
plot((Xt-Xt(1))*25, evt_lick(idx_tr)*150,'k')
plot(allTab.TongueY(1+(tr-1)*100:tr*100),'r')
% plot(allTab.TongueProba((tr-1)*100:tr*100),'r')
legend('TTL LickSensor','DLC TongueY'  )
title(['SENSOR vs DLC trial' num2str(tr)])
% ylim([0 1.5])

%% average impulsive alignement
% trtype='omi'
if trtype=='imp'
idx_imp = sort([trial.idx_errorDelay_PL_CL,...
    trial.idx_errorDelay_PR_CR,...
    trial.idx_errorDelay_PL_CR,...
    trial.idx_errorDelay_PR_CL])
elseif trtype=='cor'
    idx_imp = trial.idx_correct_R
elseif trtype=='omi'
     idx_imp = trial.idx_NoLick

end


TTL =[]; DLC = [];
plotON = 1;
for ii=2:max(size(idx_imp))
    tr=idx_imp(ii);
    idx_tr = [idx_trial_start(tr)+1:1:idx_trial_end(tr)];
    Xt =time(idx_tr); Xt(1);
    
    if plotON==1
        figure, hold on,
        plot((Xt-Xt(1))*25, evt_lick(idx_tr)*200,'k')
        plot(allTab.TongueY(((tr-1)*100)+1:tr*100),'r')
        plot(allTab.TongueX((tr-1)*100:tr*100),'b')
        legend('TTL LickSensor','DLC TongueY','DLC TongueX'  )
        title(['SENSOR vs DLC trial' num2str(tr)])
        xlabel('frames (40ms/frame)')
        ylabel('Y (pixels)')
        ylim([80 180])
        xlim([1 70])
        saveas(gcf,[SAVEFIG_folder 'LickSensor vs DLC ' trtype ' trial#' num2str(tr) '_' MouseID],'png')
    end
    
    TTL = [TTL evt_lick(idx_tr)];
    DLC = [DLC allTab.TongueY(((tr-1)*100)+1:tr*100)];
    
end

%% plot average from start trial
figure, hold on,
plot((Xt-Xt(1))*25,  mean(TTL,2)*-100,'k')
plot(mean(DLC,2)- mean(median(DLC,2)),'r')
title(['SENSOR vs DLC trial Average n=' num2str(max(size(idx_imp)))])
xlabel('frames (40ms/frame)')
ylabel('Y (pixels)')

%% plot average from 1st TTL licksensor
% plotON = 0;
% idx_imp = sort([trial.idx_errorDelay_PL_CL,...
%     trial.idx_errorDelay_PR_CR,...
%     trial.idx_errorDelay_PL_CR,...
%     trial.idx_errorDelay_PR_CL])
% for ii=2:max(size(idx_imp))
%     tr=idx_imp(ii);
%     Y_DLC = allTab.TongueY(((tr-1)*100)+1:tr*100);
%     
%     idx_tr = [idx_trial_start(tr)+1:1:idx_trial_end(tr)];
%     Xt =time(idx_tr); Xt1=Xt(1);
%     XtC(ii)= Xt(min(find(evt_lick(idx_tr)==1)));
%     XT0 = Xt -Xt1;
%     XtC0(ii) = XtC(ii)-Xt1;;
%     FrameCenter(ii) = XtC0(ii)*25;
%     idxFrameCenter(ii) = round(XtC0(ii)*25);
%     
%     XtCenter1 = (Xt-Xt1)*25;
%     XtCenter2 = (Xt-XtC(ii))*25;
%     if plotON==1
%         figure, hold on,
%         plot(XtCenter2, evt_lick(idx_tr)*200,'k');
%         plot([1:1:100]-FrameCenter(ii), Y_DLC,'r');
%     end
%     DLC(:,ii) = Y_DLC(idxFrameCenter(ii)-10:idxFrameCenter(ii)+10);
%     
% end
% 
% %%
% k=2
% figure, hold on,
% ntr = size(DLC,2)-1;
% Xframe = [-10:1:10];
% stdDLC = std(DLC)'; 
% semDLC = stdDLC(2:end)/sqrt(ntr);
% MeanDLC = mean(DLC,2) - mean(median(DLC,2));
% plotshaded(Xframe, [MeanDLC-(k*semDLC) MeanDLC  MeanDLC+(k*semDLC)],'r')
% xline(0)
% title(['SENSOR vs DLC trial Average n=' num2str(ntr)])
% xlabel('frames (40ms/frame)')
% ylabel('Y (pixels)')
