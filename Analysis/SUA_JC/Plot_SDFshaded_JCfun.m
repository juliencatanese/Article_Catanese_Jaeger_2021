function Plot_SDFshaded_JCfun(SMA, SSA, cell2plot, parfig, saveON)
% function Plot_Raster_SDFshaded_JCfun(SMA, SSA, spxtimes, trigtimes)
% param 
pre = parfig.pre;
post = parfig.post;
center_evt = parfig.center_evt;
trial_type = parfig.trial_type;
col= parfig.col;

MouseID = cell2plot.MouseID;
Day = cell2plot.Day;
ChanID = cell2plot.ChanID;
CLUST = cell2plot.CLUST;

%% PLOT RASTER

figure, 
NbTtrialType = size(SMA,1); 
for ii=1:NbTtrialType
    plot([-pre:1:post], SMA(ii,:),'color',col{ii},'LineWidth',2);
    hold on, 
    plotshaded([-pre:1:post], [SMA(ii,:)+SSA(ii,:); SMA(ii,:)-SSA(ii,:)] ,col{ii});

    pos=[-750; 0; +750]; 
    ym = max(ylim);
    YY = [0 ym];
    for ipos=1:max(size(pos))
        hold on, line([pos(ipos) pos(ipos)], YY, 'Color','k','LineStyle','--','LineWidth',1 )
    end
end
%%
ylabel(['fr (Hz)'],'FontSize',11,'FontWeight', 'normal')
xlabel(['time from ' center_evt ' start (ms)'], 'FontSize',11,'FontWeight', 'normal');

xlim([parfig.xlim])
ylim([0 ym])

title([MouseID ' ' Day ' ' ChanID ' clust#' num2str(CLUST)],'FontSize',11,'FontWeight', 'bold')

if saveON == 1
    FolderTarget = 'D:\DATA EMORY\JC_Figures\SUA\SDF_sig'
    saveas(gcf, [FolderTarget '\' MouseID '_' Day '_' ChanID '_clust#' num2str(CLUST) '_SDFshaded'],'png')
    saveas(gcf, [FolderTarget '\' MouseID '_' Day '_' ChanID '_clust#' num2str(CLUST) '_SDFshaded'],'eps')
end
