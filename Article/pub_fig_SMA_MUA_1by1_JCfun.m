function pub_fig_SMA_MUA_1by1_JCfun(idx_list, SMA, parfig)
% Fig 2A : Population firing rate
% written by Julien Catanese 3/15/2019

try
    col = parfig.col
    tit = parfig.title
catch
    col = {'k','m','g','c','y','r','b'}
    tit = 'UNDEFINED (parfig.tit)'
end


idx=idx_list;
Xt=[1:1:size(SMA,2)]-parfig.pre;

legend_all = []; idx1=[]; SMA2=[];
for ii= 1:size(idx,2);    
    figure, hold on, 
    SMA2=SMA(idx(:,ii),:);
    fstr= col{ii};
    idx1=idx(:,ii);

    if strcmp(parfig.plotshaded,'std')
        plotshaded(Xt,[nanmean(SMA2)+(nanstd(SMA2)) ;nanmean(SMA2); nanmean(SMA2)-(nanstd(SMA2));], fstr); pause(0.2)
    else
        plotshaded(Xt,[nanmean(SMA2)+(nanstd(SMA2)/sqrt(sum(idx1))) ;nanmean(SMA2); nanmean(SMA2)-(nanstd(SMA2)/sqrt(sum(idx1)));], fstr);
    end
    
    ylim(parfig.ylim)
    xlim(parfig.xlim);
    ylabel(parfig.ylabel) ;
    xlabel('time (ms)');
    title(tit);
    
    pre = 0; %parfig.pre;
%     if parfig.center_evt == 'GoCue'
%         hold on, plot([pre-1500 pre-1500], ylim,'k--','LineWidth',2);
%         hold on, plot([pre-750 pre-750], ylim,'k--', 'LineWidth',2);
%         hold on, plot([pre pre], ylim,'k--','LineWidth',2);
%     elseif parfig.center_evt == 'Delay'
%         hold on, plot([pre-750 pre-750], ylim,'k--', 'LineWidth',2);
%         hold on, plot([pre pre], ylim,'k--','LineWidth',2);
%         hold on, plot([pre+750 pre+750], ylim,'k--','LineWidth',2);
%     end
    
pause(0.1)
saveas(gcf, ['D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_FRPOP_curve' num2str(ii) ], 'png');
saveas(gcf, ['D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_FRPOP_curve' num2str(ii) ], 'emf');
pause(0.1)
    
end