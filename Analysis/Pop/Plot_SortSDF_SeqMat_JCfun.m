function [sort_sdfall] = Plot_SortSDF_SeqMat_JCfun(sdf_norm_mean_ALL, setplot, pre)

clear SNMA SNSA IDX sort_sdfall;

SNMA = sdf_norm_mean_ALL;
for irow= 1:size(SNMA,1);
    IDX(irow) =  min(find(SNMA(irow,:)==1));
end
SNMA = [IDX' SNMA];
sort_sdfall = sortrows(SNMA,'ascend');
figure, imagesc(sort_sdfall(:,2:end));
colorbar,;

colormap(setplot.colormap)
caxis(setplot.caxis);
title(setplot.title);

hold on, plot([pre-750 pre-750], ylim,'w--', 'LineWidth',2.5);
hold on, plot([pre pre], ylim,'w--','LineWidth',2.5);
hold on, plot([pre+750 pre+750], ylim,'w--','LineWidth',2.5);

ylabel ('# neurons (sorted)');
xlabel('time (ms)');

try;
    xlim(setplot.xlim);
catch;
    setplot.xlim = [pre-1500 pre+1500];
end;