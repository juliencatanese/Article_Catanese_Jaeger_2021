function Plot_SortSDF_SeqMat_JCfun(sort_mat_norm, parfig)
% function Plot_SortSDF_SeqMat_JCfun(sort_mat_norm, parfig)
% require to run "pub_comp_SortSDFMat_JCfun.m" prealably
% written by Julien Catanese 13MAr2019

figure;
imagesc(sort_mat_norm);
colorbar,;
pre=parfig.pre;

colormap(parfig.colormap);
caxis(parfig.caxis);
title(parfig.title);

if parfig.center_evt == 'GoCue';
    hold on, plot([pre-1500 pre-1500], ylim,'w--','LineWidth',2.5);
    hold on, plot([pre-750 pre-750], ylim,'w--', 'LineWidth',2.5);
    hold on, plot([pre pre], ylim,'w--','LineWidth',2.5);
elseif parfig.center_evt == 'Delay';
    hold on, plot([pre-750 pre-750], ylim,'w--', 'LineWidth',2.5);
    hold on, plot([pre pre], ylim,'w--','LineWidth',2.5);
    hold on, plot([pre+750 pre+750], ylim,'w--','LineWidth',2.5);
end

ylabel ('# neurons (sorted)');
xlabel('time (ms)');
xlim(parfig.xlim)