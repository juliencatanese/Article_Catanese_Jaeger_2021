function [VAR1norm VAR2norm VAR1sem VAR2sem ntr1 ntr2]=plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
% function plot_Abs_Average_DLCvar_JCfun(VAR1, VAR2, VARname_str1, VARname_str2)
% Julien Catanese 9/23/2020

%% Plot VAR using idxL and idxR

VAR=VAR_DLC;
TITLE_EVT = VAR_str;
Thr = mean(mean(VAR(10:30,:)));


figure, hold on, 
% initial plot order for legend 
K=1
X=[1:1:100]'
ntr1=size(idx1,2);
ntr2=size(idx2,2);
VAR1norm = nanmean(abs(VAR(:,idx1)'))-Thr;
VAR2norm = nanmean(abs(VAR(:,idx2)'))-Thr; if max(size(VAR2norm))~=100;  VAR2norm = nan(1,100); end;
VAR1sem= std(abs(VAR(:,[idx1])-Thr)',1)/sqrt(ntr1);
VAR2sem= std(abs(VAR(:,[idx2])-Thr)',1)/sqrt(ntr2); if max(size(VAR2sem))~=100;  VAR2sem = nan(1,100); end;

plot(X', VAR1norm,'r','LineWidth',2);
plot(X', VAR2norm,'b','LineWidth',2);
MeanABS = mean(abs(VAR(:,sort([idx1 idx2]))-Thr),2)
SemABS = (std(abs(VAR(:,sort([idx1 idx2]))-Thr)',1)/sqrt(ntr1+ntr2))'

plot(MeanABS,'k','LineWidth',2);
plotshaded(X', [MeanABS-(SemABS*K) MeanABS+(SemABS*K)],'k');
xline(1.5*25,'k--','LineWidth',2);
ntr1=size(idx1,2);
ntr2=size(idx2,2);
% for ii=1:ntr1; plot(VAR(:,idx1(ii)),'m'); end
% for ii=1:ntr2; plot(VAR(:,idx2(ii)),'c'); end

% replot for visual clarity
% plot(mean(VAR(:,idx1),2),'r','LineWidth',2);
% plot(mean(VAR(:,idx2),2),'b','LineWidth',2);
plot(ones(1,100)*Thr,'k--');
xline(1.5*25,'k--','LineWidth',2);

max1=max(abs(mean(VAR(:,idx1)-Thr,2))) ;
max2=max(abs(mean(VAR(:,idx2)-Thr,2)));
% ylim([Thr-max(max2,max1)-5 Thr+max(max2,max1)+5]);
ylim([-max(max2,max1)-5 max(max2,max1)+5]);

legend([idx_str1 ' (n=' num2str(ntr1) ')' ],  [idx_str2 ' (n=' num2str(ntr2) ')'],'GoCue');
title(TITLE_EVT);

% saveas(gcf, ['Average_DLCresults_' VAR_str ],'png')
% saveas(gcf, ['Average_DLCresults_' VAR_str ],'emf')
