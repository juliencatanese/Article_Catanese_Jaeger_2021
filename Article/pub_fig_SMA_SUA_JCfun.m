function pub_fig_SMA_SUA_JCfun(listc, SMA, parfig)
% function pub_fig_SMA_JCfun(listcell, parfig)
% fast plot SMA only for vizualisation of a specific list of cells
% plot SDF all 
% written by Julien Catanese the 12 March 2019

legendtext= [];
Xpos = [-parfig.pre:1:parfig.post];
figure,
for nc = 1:size(listc,1)
    listc(nc,:)
    
    hold on,
    plot(Xpos, SMA(listc.ncell(nc),:));
    
    if listc.CLUST(nc,:) <10
        legendtext= [legendtext ; [num2str(listc.MouseID(nc,:)) ' ' num2str(listc.Day(nc,:)) ' ' num2str(listc.ChanID(nc,:)) ' clu#0' num2str(listc.CLUST(nc,:))]];
    else
        legendtext= [legendtext ; [num2str(listc.MouseID(nc,:)) ' ' num2str(listc.Day(nc,:)) ' ' num2str(listc.ChanID(nc,:)) ' clu#' num2str(listc.CLUST(nc,:))]];
    end
    
end

pos=[-1500; -750; 0; +750; +1500];
YM = max(ylim);
ym = min(ylim);
YY = [ym YM];

for ipos=1:max(size(pos))
    hold on, line([pos(ipos) pos(ipos)], YY, 'Color','k','LineStyle','--','LineWidth',0.5);
end

ylabel(parfig.ylabel,'FontSize',11,'FontWeight', 'normal');
xlabel(['time from ' parfig.center_evt ' start (ms)'], 'FontSize',11,'FontWeight', 'normal');

xlim([parfig.xlim]);
% ylim([0 ym])
legend(legendtext);
title(parfig.title);

