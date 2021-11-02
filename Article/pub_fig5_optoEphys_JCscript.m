% pub_fig5_optoEphys_JCScript 
% Fig5C OPTO EPHYS RESULTS  
% by JC lastupdate 8/31/2019 
close all 
clearvars -except mypath parfig
load('listcell.mat'); 
load('Tfig2_cor.mat');
parfig.plot=0;
parfig.saveSMA=1;
parfig.saveTABLE =1;
parfig.BaselineEpoch= [150:2150]; BLE = parfig.BaselineEpoch;
parfig.pre = parfig.BaselineEpoch(end) + 1500 % define how much time before zero
parfig.post = 1500 % define how much time after zero (zero will be defined by trigtimes.cor)
parfig.center_evt =  'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)

%% compute all SMA for each trial types 
% CORRECT OPTO TRIALS
% parfig.trial_type = {'ocC'}
% subT = listcell(subT_bool,:)
% pub_comp_ephySMA_task_ttest_JCfun(subT, parfig);
% % IPSI OPTO TRIALS
% parfig.trial_type = {'oCL'}
% subT = listcell(subT_bool,:)
% pub_comp_ephySMA_task_ttest_JCfun(subT, parfig);
% % CONTRA OPTO TRIALS
% parfig.trial_type = {'oCR'}
% subT = listcell(subT_bool,:)
% pub_comp_ephySMA_task_ttest_JCfun(subT, parfig);
% %% OMISSION OPTO TRIALS
% parfig.trial_type = {'oNO'}
% subT = listcell(subT_bool,:)
% pub_comp_ephySMA_task_ttest_JCfun(subT, parfig);

%%
load('SMA_opt_GoCue360.mat');
load Tfig2_cor.mat
load Tfig3_RvL.mat
load Tfig6_opt.mat
% pub_table6_Tfig6_opt_JCscript;
% save('Tfig6_opt.mat','Tfig6_opt')

%% Define which subset of cells/sessions to use 
subT_11sess_545 = logical(Tfig2_cor.Opto_post_sess); % nSess==[2,3,4,6,7,9,10,12,13,14,15] = 360 cells
subT_11sess_360 = logical(Tfig5_imp.Opto_post_sess); % nSess==[2,3,4,6,7,9,10,12,13,14,15] = 360 cells 
subT_5sess_545 = logical(sum(Tfig2_cor.nSess == [4,7,10,13,15],2)); % subT_142cells = 5 Best Sessions with opto_trials at least 7 = 142 cells, 3 mice
subT_5sess_360 = logical(sum(Tfig5_imp.nSess == [4,7,10,13,15],2)); % subT_142cells = 5 Best Sessions with opto_trials at least 7 = 142 cells

subT_bool = subT_5sess_545; % only 142 cells are compatible with zSMA_oCL, zSMA_oCR, zSMA_oCC
subT_bool360 = subT_5sess_360; % only 142 cells are compatible with zSMA_oCL, zSMA_oCR, zSMA_oCC

%% PlotShaded Average Pop Z-score of the Firing rate during TASK : OPTO+ vs OPTO-
figure,
K=2;
load SMA_opt_GoCue360.mat  % ALL opto trials 11 sessions 
zSMA=[]; zSMA1=[]; SMA2=[]; 
SMA = SMA(subT_bool360,:); SSemA = SSemA(subT_bool360,:);
zSMA = (SMA-FRepoch.BLE.mean(subT_bool360))./mean(SSemA,2); zSMA(find(zSMA==inf))=NaN;
% zSMA = (SMA-FRepoch.BLE.mean)./FRepoch.BLE.Sem;

% to plot the thal_NOL (black curve)
idx=logical(Tfig2_cor.VMVL & Tfig2_cor.z_exct & ~Tfig2_cor.Opto_inib & ~Tfig2_cor.Opto_exct & Tfig2_cor.Opto_post_sess); 
SMA2=zSMA(idx(subT_bool360),:); col='k'; fstr= col; idx1=idx(subT_bool360);
Xt=[1:max(size(nanmean(SMA2)))]-parfig.pre;
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;nanmean(SMA2) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;  ], fstr)

% to plot the opto_inib (cyan curve)
idx=logical(Tfig2_cor.VMVL & Tfig2_cor.Opto_inib & Tfig2_cor.z_exct); 
SMA2=zSMA(idx(subT_bool360),:); col='c'; fstr= col; idx2=idx(subT_bool360);
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;nanmean(SMA2) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;  ], fstr)

% to plot the opto_exct (magenta curve)
% idx=logical(Tfig2_cor.VMVL & Tfig2_cor.Opto_exct & Tfig2_cor.z_exct); SMA2=zSMA_oCL(idx,:); col='m'; fstr= col; idx3=idx;
% hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;nanmean(SMA2) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;  ], fstr)

XoptoStim = Xt(parfig.pre-500: parfig.pre+500);
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)

ylabel('Zscore'); xlabel('time');
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);

xlim([-2500 1000]);
ylim([-2 7])
legend(['sem n=' num2str(sum(idx1)) ],' Thal NOL',['sem n=' num2str(sum(idx2))],' SvTh- ', 'Location','northwest')

%% 142   PLOT
% load('Tfig2_cor.mat');
T142=Tfig2_cor(subT_bool,:);

load SMA_oCL_GoCue142.mat
zSMA_oCL = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);  %Ntrials_oCL
load SMA_oCR_GoCue142.mat
zSMA_oCR = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);  %Ntrials_oCR
load SMA_ocC_GoCue142.mat
zSMA_oCC = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);  %Ntrials_oCC

load SMA_cCL_GoCue545.mat
SMA = SMA(subT_bool,:); SSemA = SSemA(subT_bool,:);
zSMA_cCL = (SMA-FRepoch.BLE.mean(subT_bool))./mean(SSemA,2);
load SMA_cCR_GoCue545.mat
SMA = SMA(subT_bool,:); SSemA = SSemA(subT_bool,:);
zSMA_cCR = (SMA-FRepoch.BLE.mean(subT_bool))./mean(SSemA,2);


Xt=[1:max(size(nanmean(SMA)))]-parfig.pre;
XoptoStim = Xt(parfig.pre-500: parfig.pre+500);
K=2;
idx_NOL =logical(T142.VMVL & T142.z_exct & ~T142.Opto_inib & ~T142.Opto_exct & T142.Opto_post_sess);
idx_SvTh =logical(T142.VMVL & T142.Opto_inib & T142.z_exct);

% Th- cells trials: ipsi vs opto 
figure,
SMA1=zSMA_cCL(idx_SvTh,:); col1='r'; idx1=idx_SvTh; %Ntr = Ntrials_cCL;
SMA2=zSMA_oCL(idx_SvTh,:); col2='c'; idx2=idx_SvTh; %Ntr = Ntrials_oCL;
hold on, plot(Xt, nanmean(SMA1),col1)
hold on, plot(Xt, nanmean(SMA2),col2)
hold on, plotshaded(Xt,[nanmean(SMA1)+(K*(nanstd(SMA1)/sqrt(sum(idx1)))) ; nanmean(SMA1)-(K*(nanstd(SMA2)/sqrt(sum(idx1)))) ;  ], col1)
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ;  ], col2)
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([-2500 1000]); ylim([-2 7]); ylabel('Zscore'); xlabel('time');
legend(['NONopto ipsi trials (n=' num2str(sum(idx1)) 'cells SvTh-)'],[' OPTO ipsi trials (n=' num2str(sum(idx1)) ' cells SvTh-)'], 'Location','best')
% title(['IPSI ' num2str(round(mean(Ntr)))  ' trials/sess  (+/-'  num2str(round(std(Ntr)))  ')' ] )


% Th- cells trials: contra vs opto 
figure,
SMA1=zSMA_cCR(idx_SvTh,:); col1='b'; idx1=idx_SvTh; %Ntr = Ntrials_cCR;
SMA2=zSMA_oCR(idx_SvTh,:); col2='c'; idx2=idx_SvTh; %Ntr = Ntrials_oCR;
hold on, plot(Xt, nanmean(SMA1),col1)
hold on, plot(Xt, nanmean(SMA2),col2)
hold on, plotshaded(Xt,[nanmean(SMA1)+(K*(nanstd(SMA1)/sqrt(sum(idx1)))) ; nanmean(SMA1)-(K*(nanstd(SMA2)/sqrt(sum(idx1)))) ;  ], col1)
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ;  ], col2)
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([-2500 1000]); ylim([-2 7]); ylabel('Zscore'); xlabel('time');
legend(['NONopto contra trials (n=' num2str(sum(idx1)) 'cells SvTh-)'],[' OPTO contra trials (n=' num2str(sum(idx1)) ' cells SvTh-)'], 'Location','best')
% title(['CONTRA ' num2str(round(mean(Ntr)))  ' trials/sess  (+/-'  num2str(round(std(Ntr)))  ')' ] )

%% n=142 cells PLOT Average opt trial vs non-opto 
load('Tfig2_cor.mat');
T142=Tfig2_cor(subT_bool,:);

load SMA_oNO_GoCue142.mat
zSMA_oNO = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);  
load SMA_ocC_GoCue142.mat
zSMA_oCC = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);  

load SMA_omi_GoCue360.mat
subT_bool360 = logical(sum(Tfig5_imp.nSess == [4,7,10,13,15],2));
SMA = SMA(subT_bool360,:); SSemA = SSemA(subT_bool360,:);
zSMA_omi = (SMA-FRepoch.BLE.mean(subT_bool360))./mean(SSemA,2);
load SMA_cor_GoCue545
SMA = SMA(subT_bool,:); SSemA = SSemA(subT_bool,:);
zSMA_cor = (SMA-FRepoch.BLE.mean(subT_bool))./mean(SSemA,2);

Xt=[1:max(size(nanmean(SMA)))]-parfig.pre;
XoptoStim = Xt(parfig.pre-500: parfig.pre+500);
K=2;
idx_NOL =logical(T142.VMVL & T142.z_exct & ~T142.Opto_inib & ~T142.Opto_exct & T142.Opto_post_sess);
idx_SvTh =logical(T142.VMVL & T142.Opto_inib & T142.z_exct);
idx_SvTh_pos =logical(T142.VMVL & T142.Opto_exct & T142.z_exct);

% Th- cells trials: omi vs opto   
figure,
SMA1=zSMA_omi(idx_SvTh,:); col1='g'; idx1=idx_SvTh; 
SMA2=zSMA_oNO(idx_SvTh,:); col2='c'; idx2=idx_SvTh; 
hold on, plot(Xt, nanmean(SMA1),col1)
hold on, plot(Xt, nanmean(SMA2),col2)
hold on, plotshaded(Xt,[nanmean(SMA1)+(K*(nanstd(SMA1)/sqrt(sum(idx1)))) ; nanmean(SMA1)-(K*(nanstd(SMA2)/sqrt(sum(idx1)))) ;  ], col1)
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ;  ], col2)
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([-2500 1000]); ylim([-2 7]); ylabel('Zscore'); xlabel('time');
legend(['NONopto omi trials (n=' num2str(sum(idx1)) 'cells SvTh-)'],[' OPTO omi trials (n=' num2str(sum(idx1)) ' cells SvTh-)'], 'Location','best')
% title(['IPSI ' num2str(round(mean(Ntr)))  ' trials/sess  (+/-'  num2str(round(std(Ntr)))  ')' ] )

% NOL cells trials: omi vs opto   
figure,
SMA1=zSMA_omi(idx_SvTh,:); col1='g'; idx1=idx_SvTh; 
SMA2=zSMA_oNO(idx_SvTh,:); col2='c'; idx2=idx_SvTh; 
hold on, plot(Xt, nanmean(SMA1),col1)
hold on, plot(Xt, nanmean(SMA2),col2)
hold on, plotshaded(Xt,[nanmean(SMA1)+(K*(nanstd(SMA1)/sqrt(sum(idx1)))) ; nanmean(SMA1)-(K*(nanstd(SMA2)/sqrt(sum(idx1)))) ;  ], col1)
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ;  ], col2)
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([-2500 1000]); ylim([-2 7]); ylabel('Zscore'); xlabel('time');
legend(['NONopto omi trials (n=' num2str(sum(idx1)) 'cells SvTh-)'],[' OPTO omi trials (n=' num2str(sum(idx1)) ' cells SvTh-)'], 'Location','best')
% title(['IPSI ' num2str(round(mean(Ntr)))  ' trials/sess  (+/-'  num2str(round(std(Ntr)))  ')' ] )

% Th+ cells trials: omi vs opto   
figure,
SMA1=zSMA_omi(idx_SvTh,:); col1='g'; idx1=idx_SvTh; 
SMA2=zSMA_oNO(idx_SvTh,:); col2='c'; idx2=idx_SvTh; 
hold on, plot(Xt, nanmean(SMA1),col1)
hold on, plot(Xt, nanmean(SMA2),col2)
hold on, plotshaded(Xt,[nanmean(SMA1)+(K*(nanstd(SMA1)/sqrt(sum(idx1)))) ; nanmean(SMA1)-(K*(nanstd(SMA2)/sqrt(sum(idx1)))) ;  ], col1)
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ;  ], col2)
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([-2500 1000]); ylim([-2 7]); ylabel('Zscore'); xlabel('time');
legend(['NONopto omi trials (n=' num2str(sum(idx1)) 'cells SvTh-)'],[' OPTO omi trials (n=' num2str(sum(idx1)) ' cells SvTh-)'], 'Location','best')
% title(['IPSI ' num2str(round(mean(Ntr)))  ' trials/sess  (+/-'  num2str(round(std(Ntr)))  ')' ] )

% Th- cells trials: cor vs opto    
figure,
SMA1=zSMA_cor(idx_SvTh,:); col1='k'; idx1=idx_SvTh; 
SMA2=zSMA_oCC(idx_SvTh,:); col2='c'; idx2=idx_SvTh; 
hold on, plot(Xt, nanmean(SMA1),col1)
hold on, plot(Xt, nanmean(SMA2),col2)
hold on, plotshaded(Xt,[nanmean(SMA1)+(K*(nanstd(SMA1)/sqrt(sum(idx1)))) ; nanmean(SMA1)-(K*(nanstd(SMA2)/sqrt(sum(idx1)))) ;  ], col1)
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ;  ], col2)
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([-2500 1000]); ylim([-2 7]); ylabel('Zscore'); xlabel('time');
legend(['NONopto cor trials (n=' num2str(sum(idx1)) 'cells SvTh-)'],[' OPTO cor trials (n=' num2str(sum(idx1)) ' cells SvTh-)'], 'Location','best')
% title(['CONTRA ' num2str(round(mean(Ntr)))  ' trials/sess  (+/-'  num2str(round(std(Ntr)))  ')' ] )

% NOL cells trials: cor vs opto    
figure,
SMA1=zSMA_cor(idx_NOL,:); col1='k'; idx1=idx_NOL; 
SMA2=zSMA_oCC(idx_NOL,:); col2='c'; idx2=idx_NOL; 
hold on, plot(Xt, nanmean(SMA1),col1)
hold on, plot(Xt, nanmean(SMA2),col2)
hold on, plotshaded(Xt,[nanmean(SMA1)+(K*(nanstd(SMA1)/sqrt(sum(idx1)))) ; nanmean(SMA1)-(K*(nanstd(SMA2)/sqrt(sum(idx1)))) ;  ], col1)
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ;  ], col2)
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([-2500 1000]); ylim([-2 7]); ylabel('Zscore'); xlabel('time');
legend(['NONopto cor trials (n=' num2str(sum(idx1)) 'cells NOL)'],[' OPTO cor trials (n=' num2str(sum(idx1)) ' cells NOL)'], 'Location','best')
% title(['CONTRA ' num2str(round(mean(Ntr)))  ' trials/sess  (+/-'  num2str(round(std(Ntr)))  ')' ] )

% NOL cells trials: cor vs opto    
figure,
SMA1=zSMA_cor(idx_SvTh_pos,:); col1='k'; idx1=idx_SvTh_pos; 
SMA2=zSMA_oCC(idx_SvTh_pos,:); col2='c'; idx2=idx_SvTh_pos; 
hold on, plot(Xt, nanmean(SMA1),col1)
hold on, plot(Xt, nanmean(SMA2),col2)
hold on, plotshaded(Xt,[nanmean(SMA1)+(K*(nanstd(SMA1)/sqrt(sum(idx1)))) ; nanmean(SMA1)-(K*(nanstd(SMA2)/sqrt(sum(idx1)))) ;  ], col1)
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ;  ], col2)
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([-2500 1000]); ylim([-2 7]); ylabel('Zscore'); xlabel('time');
legend(['NONopto cor trials (n=' num2str(sum(idx1)) 'cells SvTh+)'],[' OPTO cor trials (n=' num2str(sum(idx1)) ' cells SvTh+)'], 'Location','best')
% title(['CONTRA ' num2str(round(mean(Ntr)))  ' trials/sess  (+/-'  num2str(round(std(Ntr)))  ')' ] )

%%
idx_NOL =logical(T142.VMVL & T142.z_exct & ~T142.Opto_inib & ~T142.Opto_exct & T142.Opto_post_sess);
idx_SvTh =logical(T142.VMVL & T142.Opto_inib & T142.z_exct);

T142=Tfig3_RvL(subT_bool,:); 
idx_SvTh_ipsi =logical(T142.VMVL & T142.dzRL_inib & T142.Opto_inib);
idx_SvTh_cont =logical(T142.VMVL & T142.dzRL_exct & T142.Opto_inib);

figure,

% subplot(2,2,1)
idx1=idx_SvTh_ipsi; idx2=idx_SvTh_ipsi;
SMA1=zSMA_cCL(idx1,:); col1='r'; %Ntr = Ntrials_cCL;
SMA2=zSMA_oCL(idx2,:); col2='c'; %Ntr = Ntrials_oCL;
hold on, plot(Xt, nanmean(SMA1),col1)
hold on, plot(Xt, nanmean(SMA2),col2)
hold on, plotshaded(Xt,[nanmean(SMA1)+(K*(nanstd(SMA1)/sqrt(sum(idx1)))) ; nanmean(SMA1)-(K*(nanstd(SMA2)/sqrt(sum(idx1)))) ;  ], col1)
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ;  ], col2)
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([-2500 1000]); ylim([-2 7]); ylabel('Zscore'); xlabel('time');
legend('NONopto ipsi tr' , 'OPTO ipsi tr'  , [ 'n= ' num2str(sum(idx1)) ' ipsi SvTh-' ], 'Location','best');
% title(['IPSI ' num2str(round(mean(Ntr)))  ' trials/sess  (+/-'  num2str(round(std(Ntr)))  ')' ] )

figure,
% subplot(2,2,2)
idx1=idx_SvTh_ipsi; idx2=idx_SvTh_ipsi;
SMA1=zSMA_cCR(idx1,:); col1='b'; %Ntr = Ntrials_cCR;
SMA2=zSMA_oCR(idx2,:); col2='c'; %Ntr = Ntrials_oCR;
hold on, plot(Xt, nanmean(SMA1),col1)
hold on, plot(Xt, nanmean(SMA2),col2)
hold on, plotshaded(Xt,[nanmean(SMA1)+(K*(nanstd(SMA1)/sqrt(sum(idx1)))) ; nanmean(SMA1)-(K*(nanstd(SMA2)/sqrt(sum(idx1)))) ;  ], col1)
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ;  ], col2)
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([-2500 1000]); ylim([-2 7]); ylabel('Zscore'); xlabel('time');
legend('NONopto contra tr'  , 'OPTO contra tr'  , ['n= ' num2str(sum(idx1)) ' ipsi SvTh-'], 'Location','best');

figure,
% subplot(2,2,3)
idx1=idx_SvTh_cont; idx1(46)=0; idx2=idx1;
SMA1=zSMA_cCL(idx1,:); col1='r'; %Ntr = Ntrials_cCR;
SMA2=zSMA_oCL(idx2,:); col2='c'; %Ntr = Ntrials_oCR;
hold on, plot(Xt, nanmean(SMA1),col1)
hold on, plot(Xt, nanmean(SMA2),col2)
hold on, plotshaded(Xt,[nanmean(SMA1)+(K*(nanstd(SMA1)/sqrt(sum(idx1)))) ; nanmean(SMA1)-(K*(nanstd(SMA2)/sqrt(sum(idx1)))) ;  ], col1)
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ;  ], col2)
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([-2500 1000]); ylim([-2 7]); ylabel('Zscore'); xlabel('time');
legend('NONopto ipsi tr'  , 'OPTO ipsi tr'  , ['n= ' num2str(sum(idx1)) ' contra SvTh-'], 'Location','best');
% title(['CONTRA ' num2str(round(mean(Ntr)))  ' trials/sess  (+/-'  num2str(round(std(Ntr)))  ')' ] )

figure,
% subplot(2,2,4)
idx1=idx_SvTh_cont; idx1(46)=0; idx2=idx1;
SMA1=zSMA_cCR(idx1,:); col1='b'; %Ntr = Ntrials_cCR;
SMA2=zSMA_oCR(idx2,:); col2='c'; %Ntr = Ntrials_oCR;
hold on, plot(Xt, nanmean(SMA1),col1)
hold on, plot(Xt, nanmean(SMA2),col2)
hold on, plotshaded(Xt,[nanmean(SMA1)+(K*(nanstd(SMA1)/sqrt(sum(idx1)))) ; nanmean(SMA1)-(K*(nanstd(SMA2)/sqrt(sum(idx1)))) ;  ], col1)
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx2)))) ;  ], col2)
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([-2500 1000]); ylim([-2 7]); ylabel('Zscore'); xlabel('time');
legend('NONopto contra tr'  , 'OPTO contra tr'  , ['n= ' num2str(sum(idx1)) ' contra SvTh-'], 'Location','best')
% title(['CONTRA ' num2str(round(mean(Ntr)))  ' trials/sess  (+/-'  num2str(round(std(Ntr)))  ')' ] )

%%

figure,
Xopt = parfig.pre-500:parfig.pre+500

subplot(2,2,1)
data = [ mean(zSMA_cCL(idx_SvTh_ipsi,Xopt)-zSMA_oCL(idx_SvTh_ipsi,Xopt),2) ]
bar(data); ylim([-3 3]);
subplot(2,2,2)
data = [ mean(zSMA_cCR(idx_SvTh_ipsi,Xopt)-zSMA_oCR(idx_SvTh_ipsi,Xopt),2) ]
bar(data); ylim([-3 3]);
subplot(2,2,3)
data = [ mean(zSMA_cCL(idx_SvTh_cont,Xopt)-zSMA_oCL(idx_SvTh_cont,Xopt),2) ]
bar(data); ylim([-3 3]);
subplot(2,2,4)
data = [ mean(zSMA_cCR(idx_SvTh_cont,Xopt)-zSMA_oCR(idx_SvTh_cont,Xopt),2) ]% ;  ;  mean(SMA1(:,Xopt)-SMA2(:,Xopt)) ;  mean(SMA1(:,Xopt)-SMA2(:,Xopt));  mean(SMA1(:,Xopt)-SMA2(:,Xopt))]
bar(data); ylim([-3 3]);

figure,
subplot(2,2,1)
data = [ mean(zSMA_cCL(idx_SvTh_ipsi,Xopt)-zSMA_oCL(idx_SvTh_ipsi,Xopt)) ]
bar(data); ylim([-3 3]);
subplot(2,2,2)
data = [ mean(zSMA_cCR(idx_SvTh_ipsi,Xopt)-zSMA_oCR(idx_SvTh_ipsi,Xopt)) ]
bar(data); ylim([-3 3]);
subplot(2,2,3)
data = [ mean(zSMA_cCL(idx_SvTh_cont,Xopt)-zSMA_oCL(idx_SvTh_cont,Xopt)) ]
bar(data); ylim([-3 3]);
subplot(2,2,4)
data = [ mean(zSMA_cCR(idx_SvTh_cont,Xopt)-zSMA_oCR(idx_SvTh_cont,Xopt)) ]% ;  ;  mean(SMA1(:,Xopt)-SMA2(:,Xopt)) ;  mean(SMA1(:,Xopt)-SMA2(:,Xopt));  mean(SMA1(:,Xopt)-SMA2(:,Xopt))]
bar(data); ylim([-3 3]);

%%
figure,
data = [mean(mean(zSMA_cCL(idx_SvTh_ipsi,Xopt)-zSMA_oCL(idx_SvTh_ipsi,Xopt),2)),...
    mean(mean(zSMA_cCR(idx_SvTh_ipsi,Xopt)-zSMA_oCR(idx_SvTh_ipsi,Xopt),2)),...
    mean(mean(zSMA_cCL(idx_SvTh_cont,Xopt)-zSMA_oCL(idx_SvTh_cont,Xopt),2)),...
    mean(mean(zSMA_cCR(idx_SvTh_cont,Xopt)-zSMA_oCR(idx_SvTh_cont,Xopt),2)) ]
bar(data); ylim([-3 3]);
disp('ttest2 2 tails')
[h p ] = ttest2(mean(zSMA_cCL(idx_SvTh_ipsi,Xopt),2), mean(zSMA_oCL(idx_SvTh_ipsi,Xopt),2)); % p=0.0057
[h p ] = ttest2(mean(zSMA_cCR(idx_SvTh_ipsi,Xopt),2), mean(zSMA_oCR(idx_SvTh_ipsi,Xopt),2))
[h p ] = ttest2(mean(zSMA_cCL(idx_SvTh_cont,Xopt),2), mean(zSMA_oCL(idx_SvTh_cont,Xopt),2))
[h p ] = ttest2(mean(zSMA_cCR(idx_SvTh_cont,Xopt),2), mean(zSMA_oCR(idx_SvTh_cont,Xopt),2)); % p=0.0341

disp('ranksum 2 tails ')
[ p h] = ranksum(mean(zSMA_cCL(idx_SvTh_ipsi,Xopt),2), mean(zSMA_oCL(idx_SvTh_ipsi,Xopt),2)) % p=0.0057
[ p h] = ranksum(mean(zSMA_cCR(idx_SvTh_ipsi,Xopt),2), mean(zSMA_oCR(idx_SvTh_ipsi,Xopt),2))
[ p h] = ranksum(mean(zSMA_cCL(idx_SvTh_cont,Xopt),2), mean(zSMA_oCL(idx_SvTh_cont,Xopt),2))
[ p h] = ranksum(mean(zSMA_cCR(idx_SvTh_cont,Xopt),2), mean(zSMA_oCR(idx_SvTh_cont,Xopt),2))% p=0.0569

disp('ranksum 1 tail ')
[ p h] = ranksum(mean(zSMA_cCL(idx_SvTh_ipsi,Xopt),2), mean(zSMA_oCL(idx_SvTh_ipsi,Xopt),2),'tail','right') % p=0.0043
[ p h] = ranksum(mean(zSMA_cCR(idx_SvTh_ipsi,Xopt),2), mean(zSMA_oCR(idx_SvTh_ipsi,Xopt),2),'tail','right')
[ p h] = ranksum(mean(zSMA_cCL(idx_SvTh_cont,Xopt),2), mean(zSMA_oCL(idx_SvTh_cont,Xopt),2),'tail','right')
[ p h] = ranksum(mean(zSMA_cCR(idx_SvTh_cont,Xopt),2), mean(zSMA_oCR(idx_SvTh_cont,Xopt),2),'tail','right')% p=0.0284

