% Perc_STATpeak_distrib_ALLSession_JCscript
% need to use SDF_mean_sem_Baseline_mlibJCscript
% Written by Julien Catanese 10/27/2018

%% Define Session list
clear all
close all    

figplot=0;

cd('D:\JC_Analysis');
SessList = dir(['*/*taskopto*']);

NSess= max(size(SessList)); % Number of Sessions
ALL_SDFmean = [];ALL_SDFsem = []; cCR=[]; cCL=[]; 

%% loop trhough all Sessions named "taskopto"
for ns=1:NSess
    SessID= SessList(ns).name;
    SessPath=SessList(ns).folder;
    cd([SessPath '\' SessID]);
    
    % Set Parameters
    psth_trig_evt =  'Delay'; % center for SDF
    psth_trial_type = {'cor'};
    pre=750+2250; % 2250ms baseline + 750ms puff
    post=1500; %750ms delay + 750ms response

    SDF_mean_sem_Baseline_mlibJCscript
    
    % create the ALL-Matrice
    ALL_SDFmean = [ALL_SDFmean; sort_sdfall] ;
    ALL_SDFsem = [ALL_SDFsem; sort_sdfall_sem];
    
end
%%
ALL_SDFmean_sorted = sortrows(ALL_SDFmean,'ascend');
ALL_SDFmean_sort = ALL_SDFmean_sorted(:,2:end);

ALL_SDFsem_sorted = sortrows(ALL_SDFsem,'ascend');
ALL_SDFsem_sort = ALL_SDFsem_sorted(:,2:end);

figure, imagesc(ALL_SDFmean_sort);
colorbar,;
caxis([0.5 1]);

%% Compute Average and SEM of the Baseline firing rate

ncell=size(ALL_SDFmean_sort,1);
AllB = [];
close all
for i=1:ncell
    
    Base1 = ALL_SDFmean_sort(i,163:2250); %  remove the first 163ms for edge effect (conv Gauss Kernel)
    Base2 = ALL_SDFsem_sort(i,163:2250);
    Base1_Av= mean(Base1)
    Base2_sem = mean(Base2)
    
    thr= Base1_Av + (3*Base2_sem)
    Cell = ALL_SDFmean_sort(i,:)

    if figplot==1
        figure(i)
        fstr='m';
        X = [-pre:1:post];
        Y = ones(1,pre+post+1) ;
        hold on, plot(X, Cell );
        hold on, plot(X, Y*Base1_Av, 'm')
        hold on, plotshaded(X,[Y*(Base1_Av+Base2_sem); Y*(Base1_Av-Base2_sem)],fstr)
        hold on, plot(X, Y*thr, 'k--');
        
        legend('SDF','Baseline (mean fr)','sem','trh (Baseline+3sem)');
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
figure,
ncell=i
X=[-pre:1:post]
Y=sum(AllB)/ncell
w=25
Gauss_width = max([11 6*w+1]); % hm, should be an odd number... e.g. 11
kernel      = normpdf(-floor(Gauss_width/2):floor(Gauss_width/2),0,w);
hold on, plot(X, conv(Y,kernel,'same'),'Color',[0.3 0.4 0.5],'LineWidth',2.5)
% hold on, plot(X, Y,'r','LineWidth',0.5)
ylabel('% modulation');
xlabel('time from delay start (ms)');
title (['distribution of modulation over #' num2str(ncell) ' cells'])

xlim([-1500 1500])
ylim([0.0 0.7])
hold on, plot([-750 -750], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5)
hold on, plot([0 0], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5)
hold on, plot([750 750], ylim,'Color','k','LineStyle','--' ,'LineWidth',1.5)
