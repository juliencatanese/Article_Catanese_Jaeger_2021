function [Yconv, Y] = pub_fig_distrib_percZthr_JCfun(BarSig, parfig)
% function plot_distrib_Peak3sem_JCfun(BarSig, ncell, col, pre, post)
% Make a plot of the distribution of the percentage of modulated activity in time
% parfig.title = ['Distribution (#Cell='  num2str(nsc) ' Zthr >' num2str(Zthr) ' trials-' parfig.trial_type{1} ];
% parfig.xlim = [0-2250 0+750];
% parfig.ylim = [0  50]
% parfig.Zthr = 10;
% writen by julien Catanese 12/05/2018

ncell=size(BarSig,1);
% define Axis
X=[-parfig.pre:1:parfig.post];
Y=(sum(BarSig)/ncell)*100;

% create a Gaussian Kernel
w=75; %ms
Gauss_width = max([11 6*w+1]); % hm, should be an odd number... e.g. 11
kernel      = normpdf(-floor(Gauss_width/2):floor(Gauss_width/2),0,w);
Yconv =  conv(Y,kernel,'same');


% Plot the smoothed histogram
if parfig.plot == 1
    hold on, plot(X, Yconv,'LineWidth',2.5);
    title(parfig.title);
    xlim(parfig.xlim);
    ylim(parfig.ylim);
    ylabel('% cell modulated Above thr');
    xlabel('time from GoCue (ms)');
end

% dashed lines to separates epochs
if parfig.plotEpochLines == 1
    pre=0;
    if parfig.center_evt == 'GoCue'
        hold on, plot([pre-1500 pre-1500], ylim,'k--','LineWidth',2.5);
        hold on, plot([pre-750 pre-750], ylim,'k--', 'LineWidth',2.5);
        hold on, plot([pre pre], ylim,'k--','LineWidth',2.5);
    elseif parfig.center_evt == 'Delay'
        hold on, plot([pre-750 pre-750], ylim,'k--', 'LineWidth',2.5);
        hold on, plot([pre pre], ylim,'k--','LineWidth',2.5);
        hold on, plot([pre+750 pre+750], ylim,'k--','LineWidth',2.5);
    end
end
