function plot_distrib_3zscore_JCfun(BarSig, parfig)
% function plot_distrib_Peak3sem_JCfun(BarSig, ncell, col, pre, post)
% Make a plot of the distribution of the percentage of modulated activity in time 
% writen by julien Catanese 12/05/2018

pre= parfig.pre;
post=parfig.post
% XLim= parfig.xlim;
TITLE = parfig.title;
col=parfig.col;

% define Axis 
X=[-pre:1:post];
Y=(sum(BarSig)/ncell)*100;

% create a Gaussian Kernel 
w=25%ms
Gauss_width = max([11 6*w+1]); % hm, should be an odd number... e.g. 11
kernel      = normpdf(-floor(Gauss_width/2):floor(Gauss_width/2),0,w);

% Plot the smoothed histogram 
hold on, plot(X, conv(Y,kernel,'same'),'Color',col{1} ,'LineWidth',2.5);

% labels 
ylabel('% modulation');
xlabel('time from delay start (ms)');
title(TITLE);
% xlim([XLim]);

% dashed lines to separates epochs 
hold on, plot([-750 -750], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5);
hold on, plot([0 0], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5);
hold on, plot([750 750], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5);


