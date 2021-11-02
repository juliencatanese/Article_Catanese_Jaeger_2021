function Plot_SortSDF_SeqMat_JCfun(sort_mat_norm, setplot)
% function Plot_SortSDF_SeqMat_JCfun(sort_mat_norm, setplot)
% need to use 
% written by Julien Catanese 13MAr2019

figure, imagesc(sort_mat_norm);
colorbar,;
pre=setplot.pre;

colormap(setplot.colormap)
caxis(setplot.caxis);
title(setplot.title);

hold on, plot([pre-750 pre-750], ylim,'w--', 'LineWidth',2.5);
hold on, plot([pre pre], ylim,'w--','LineWidth',2.5);
hold on, plot([pre+750 pre+750], ylim,'w--','LineWidth',2.5);

ylabel ('# neurons (sorted)');
xlabel('time (ms)');
