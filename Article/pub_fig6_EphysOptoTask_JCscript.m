% pub_fig_EphysOptoTask_JCscript

%% Fig6: PLot OPTO TASK effect:  opto trials SMA vs Corr trials SMA   
% parfig.BaselineEpoch= [150:2150]; BLE = parfig.BaselineEpoch;
% parfig.center_evt =  'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
% parfig.trial_type = {'opt'};
% parfig.col = {'k'};
% parfig.xlim = [-1500 1500];
% parfig.plot=0;
% 
% subT= [];
% fig
% % % listcell2= listcell(Tcombo.Opto_inib & Tcombo.VMVL,:);
% [SMA FRepoch ttestEpoch trigtimes spxtimes sdf_alltr] = pub_comp_ephySMA_task_ttest_JCfun(subT, parfig);
% 

%% PlotShaded Average Pop Z-score of the Firing rate during TASK : OPTO+ vs OPTO-
% close all
figure,
K=2; 
load SMA_opt_GoCue545.mat % to compare different cells types during opto trial
load Tcombo
load Tfig6_opt.mat 

zSMA=[]; SMA1=[]; SMA2=[];
zSMA = (SMA-FRepoch.BLE.mean)./FRepoch.BLE.Sem;
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2); zSMA(find(zSMA==inf))=NaN;
SMA1=zSMA;

idx=logical(Tcombo.VMVL & Tfig6_opt.z_exct & ~Tcombo.Opto_inib & ~Tcombo.Opto_exct & Tcombo.Opto_post_sess); SMA2=SMA1(idx,:); col='k'; fstr= col; idx1=idx;
X=[1:max(size(nanmean(SMA2)))]-parfig.pre;
hold on, plotshaded(X,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;nanmean(SMA2) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;  ], fstr)

idx=logical(Tcombo.VMVL & Tcombo.Opto_inib & Tcombo.z_exct); SMA2=SMA1(idx,:); col='c'; fstr= col; idx2=idx;
hold on, plotshaded(X,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;nanmean(SMA2) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;  ], fstr)

idx=logical(Tcombo.VMVL & Tcombo.Opto_exct & Tcombo.z_exct); SMA2=SMA1(idx,:); col='m'; fstr= col; idx3=idx;
hold on, plotshaded(X,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;nanmean(SMA2) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;  ], fstr)

XoptoStim = X(parfig.pre-500: parfig.pre+500);
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)

ylabel('Zscore'); xlabel('time');
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);

xlim([-2500 1000]); 
ylim([-1 5])
legend(['sem n=' num2str(sum(idx1)) ],' Thal NOL',['sem n=' num2str(sum(idx2))],' SvTh- ',['sem n=' num2str(sum(idx3))],'SvTh+ ', 'Location','northwest')


%% PlotShaded Average Pop Z-score of the Firing rate during TASK : IPSI vs CONTRA
close all
figure,

zSMA=[]; SMA1=[];
zSMA = (SMA-FRepoch.BLE.mean)./FRepoch.BLE.Sem;
zSMA(find(zSMA==inf))=NaN;
SMA1=zSMA;
X=(1:max(size(nanmean(SMA1))))-parfig.pre;

subplot(2,1,1)

col='b' ;SMA2=[];
idx=logical(Tcombo.VMVL & Tcombo.dzRL_exct & (~Tcombo.Opto_inib) & (~Tcombo.Opto_exct)); SMA2=SMA1(idx,:); fstr= col; idx1=idx;
hold on, plotshaded(X,[nanmean(SMA2)+(nanstd(SMA2)/sqrt(sum(idx))) ;nanmean(SMA2) ; nanmean(SMA2)-(nanstd(SMA2)/sqrt(sum(idx))) ;  ], fstr)

col='r' ;SMA2=[];
idx=logical(Tcombo.VMVL & Tcombo.dzRL_inib & ~Tcombo.Opto_inib & ~Tcombo.Opto_exct); SMA2=SMA1(idx,:); fstr= col; idx2=idx;
hold on, plotshaded(X,[nanmean(SMA2)+(nanstd(SMA2)/sqrt(sum(idx))) ;nanmean(SMA2) ; nanmean(SMA2)-(nanstd(SMA2)/sqrt(sum(idx))) ;  ], fstr)

legend( ['sem n=' num2str(sum(idx1))],'(non opto) CONTRA cells ',...
    ['sem n=' num2str(sum(idx2))],'(non opto) IPSI cells ', 'Location','northwest')
xlim([500 4500]); ylim([-3 12])
ylabel('Zscore'); xlabel('time');

subplot(2,1,2)

% Delta CONTRA - Control
col='b' ;SMA2=[];
idx=logical(Tcombo.VMVL &  Tcombo.dzRL_exct & (Tcombo.Opto_inib)); SMA2=SMA1(idx,:);  fstr= col; idx3=idx;
SMActr=SMA1(idx1,:);
DelatSMA1 = nanmean(SMA2)-nanmean(SMActr)
hold on, plot(X, DelatSMA1, 'color', col)
% Delta IPSI - Control
col='r' ;SMA2=[];
idx=logical(Tcombo.VMVL &  Tcombo.dzRL_inib & (Tcombo.Opto_inib)); SMA2=SMA1(idx,:); fstr= col; idx4=idx;
SMActr=SMA1(idx2,:);
DelatSMA2 = nanmean(SMA2)-nanmean(SMActr)
hold on, plot(X, DelatSMA2, 'color', col)

% Delta DeltaIPSI - DeltaCONTRA
col='g' ;SMA2=[];
hold on, plot(X, DelatSMA2-DelatSMA1, 'color', col, 'LineWidth' , 3)

hold on, plot(X, zeros(size(DelatSMA2)), '--k')
XoptoStim = [-500:500]
hold on, plot(XoptoStim , zeros(size(DelatSMA2(XoptoStim)))-9, 'c','LineWidth' , 3)

legend(['CONTRA Delta = (OptoContra - NonOptoContra)  n=' num2str(sum(idx3)) ' OPTOcontra cells'],...
    ['  IPSI Delta = (OptoIpsi   - NonOptoIpsi  )  n=' num2str(sum(idx4)) ' OPTOipsi cells  '],...
    [' Delta  (IPSIdelta - CONTRAdelta) '],...
    'Location','northwest')
xlim([500 4500]); ylim([-10 10])
ylabel('Delta Z Ipsi-Contra'); xlabel('time');

%% plot shaded for 69 opto inib cells with opto response in Opto trials vs corr trials
load('listcell.mat');

col='k'; SMA2=[];
load('SMA_cor_GoCue.mat');
idx=logical(Tcombo.VMVL & (Tcombo.RespCell | Tcombo.BothCell) & Tcombo.ipsi_cell);
SMA2=SMA(idx,:); fstr=col; idx1=idx; ntr_cor = size(SMA2,1)
SMAcor=nanmean(SMA2);  SSemAcor = nanstd(SMA2)/sqrt(size(SMA2,1));
X=1:max(size(nanmean(SMA2))); Xms=X-parfig.pre;
figure, hold on,  plotshaded(Xms, [SMAcor-SSemAcor; SMAcor ;SMAcor+SSemAcor],  fstr),

col='c'; SMA2=[];
load('SMA_opt_GoCue.mat'); idx=[];
idx=logical(Tcombo.VMVL & (Tcombo.RespCell | Tcombo.BothCell) & Tcombo.ipsi_cell);
SMA2=SMA(idx,:); fstr=col; idx2=idx; ntr_opt = size(SMA2,1);
SMAopt=nanmean(SMA2); SSemAopt = nanstd(SMA2)/sqrt(size(SMA2,1));
Xms=X-parfig.pre;
hold on,  plotshaded(Xms, [SMAopt-SSemAopt; SMAopt ;SMAopt+SSemAopt],  fstr),

XoptoStim = [-500:+500]
hold on, plot(XoptoStim , zeros(size((XoptoStim)))+10, 'c','LineWidth' , 3)
xlim([-2750 1000])
% ylim([5 11])
legend( ['sem n=' num2str(sum(idx1))],'Corr Trials (opto inib cells) ',...
    ['sem n=' num2str(sum(idx2))],'Opto Trials (opto inib cells)  ','Laser Stim', 'Location','northwest')
xlabel('time (sample)'), ylabel('firing rate')
