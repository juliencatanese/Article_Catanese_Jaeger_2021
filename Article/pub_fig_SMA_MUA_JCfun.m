function pub_fig_SMA_MUA_JCfun(idx_list, SMA, parfig)
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
figure,
Xt=[1:1:size(SMA,2)]-parfig.pre;

% LegText={'ALL VMVL', 'Exc VMVL' , 'Inh VMVL'}
legend_all = []; idx1=[]; SMA2=[];
for ii= 1:size(idx,2);
    SMA2=SMA(idx(:,ii),:);
    fstr= col{ii};
    idx1=idx(:,ii);
    hold on,
    k= parfig.k; 
    
    if strcmp(parfig.plotshaded,'std')
        if ii~=1
            plotshaded(Xt,[nanmean(SMA2)+(nanstd(SMA2)/sqrt(sum(idx1))) ;nanmean(SMA2); nanmean(SMA2)-(nanstd(SMA2)/sqrt(sum(idx1)));], fstr);
        elseif ii==1
            plotshaded(Xt,[nanmean(SMA2)+(nanstd(SMA2)) ;nanmean(SMA2); nanmean(SMA2)-(nanstd(SMA2));], fstr);
        end
    else
        plotshaded(Xt,[nanmean(SMA2)+(k*(nanstd(SMA2)/sqrt(sum(idx1)))) ;nanmean(SMA2); nanmean(SMA2)-(k*(nanstd(SMA2)/sqrt(sum(idx1))));], fstr);
    end
end

xlim(parfig.xlim);
ylabel(parfig.ylabel) ;
xlabel('time (ms)');
title(tit);

% pre = 0; %parfig.pre;
% if parfig.center_evt == 'GoCue'
%     hold on, plot([pre-1500 pre-1500], ylim,'k--','LineWidth',2.5);
%     hold on, plot([pre-750 pre-750], ylim,'k--', 'LineWidth',2.5);
%     hold on, plot([pre pre], ylim,'k--','LineWidth',2.5);
% elseif parfig.center_evt == 'Delay'
%     hold on, plot([pre-750 pre-750], ylim,'k--', 'LineWidth',2.5);
%     hold on, plot([pre pre], ylim,'k--','LineWidth',2.5);
%     hold on, plot([pre+750 pre+750], ylim,'k--','LineWidth',2.5);
% end

% legend('All VMVL', '' , 'Exc VMVL' , '' , 'Inh VMVL' , '' )
