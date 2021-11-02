function plot_Average_trialtype_DLCvar(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2)
% function plot_Average_trialtype_DLCvar(VAR1, VAR2, VARname_str1, VARname_str2)
% Julien Catanese 9/23/2020

%% Plot VAR using idxL and idxR

VAR=VAR_DLC;
TITLE_EVT = VAR_str;
Thr = median(mean(VAR));


figure, hold on, 
% initial plot order for legend 
plot(mean(VAR(:,idx1),2),'r','LineWidth',2);
plot(mean(VAR(:,idx2),2),'b','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
ntr1=size(idx1,2);
ntr2=size(idx2,2);
for ii=1:ntr1; plot(VAR(:,idx1(ii)),'m'); end
for ii=1:ntr2; plot(VAR(:,idx2(ii)),'c'); end

% replot for visual clarity
plot(mean(VAR(:,idx1),2),'r','LineWidth',2);
plot(mean(VAR(:,idx2),2),'b','LineWidth',2);
plot(ones(1,100)*Thr,'k--');
xline(1.5*25,'k--','LineWidth',2);

max1=max(abs(mean(VAR(:,idx1)-Thr,2))) ;
max2=max(abs(mean(VAR(:,idx2)-Thr,2)));
ylim([Thr-max(max2,max1)-5 Thr+max(max2,max1)+5]);

legend([idx_str1 ' (n=' num2str(ntr1) ')' ],  [idx_str2 ' (n=' num2str(ntr2) ')'],'GoCue');
title(TITLE_EVT);

