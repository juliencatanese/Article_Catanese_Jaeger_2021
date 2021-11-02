function plot_distrib_Peak3sem_JCfun(SMA, SVA, K, stype, col, pre, post, psth_trial_type, figplot)

%% PLot Distribution
% figure, 
ncell=i
X=[-pre:1:post]
Y=(sum(AllB)/ncell)*100
w=25
Gauss_width = max([11 6*w+1]); % hm, should be an odd number... e.g. 11
kernel      = normpdf(-floor(Gauss_width/2):floor(Gauss_width/2),0,w);
hold on, plot(X, conv(Y,kernel,'same'),'Color',col ,'LineWidth',2.5)
% hold on, plot(X, Y,'r','LineWidth',0.5)
ylabel('% modulation');
xlabel('time from delay start (ms)');
title (['distribution of peak above ' num2str(K) stype ' over #' num2str(ncell) ' cells for ' psth_trial_type{1} ' trials'])

xlim([-1500 1500])
% ylim([0 70])
hold on, plot([-750 -750], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5)
hold on, plot([0 0], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5)
hold on, plot([750 750], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5)