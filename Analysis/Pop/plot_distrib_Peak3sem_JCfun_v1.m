function plot_distrib_Peak3sem_JCfun(SMA, SVA, K, stype, col, pre, post, psth_trial_type, figplot)

ncell=size(SMA,1)
AllB = [];
for i=1:ncell
    
    Base_mean = SMA(i,163:2250); %  remove the first 163ms for edge effect (conv Gauss Kernel)
    Base_var = SVA(i,163:2250);
    
    Base_mean_Av= mean(Base_mean);
    Base_var_Av = mean(Base_var);
    
    thr= Base_mean_Av + (K*Base_var_Av);

    Cell = SMA(i,:);
    
    if figplot==1
        figure(i);
        fstr='m';
        X = [-pre:1:post];
        Y = ones(1,pre+post+1) ;
        hold on, plot(X, Cell );
        hold on, plot(X, Y*Base_mean_Av, 'm');
        hold on, plotshaded(X,[Y*(Base_mean_Av+Base_var_Av); Y*(Base_mean_Av-Base_var_Av)],fstr);
        hold on, plot(X, Y*thr, 'k--');
        
        legend('SDF','Baseline (mean fr)', stype ,['trh (Baseline+' num2str(K) stype  ')']);
        ylabel('fr norm');
        xlabel('time (ms)');
        title(['cell #' num2str(i)]);
        pause(1);
    end
    
    idxB = find(Cell>=thr);
    BarSig= zeros(1,pre+post+1);
    BarSig(idxB)= 1;
    
    AllB = [AllB; BarSig];
    
    
end

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