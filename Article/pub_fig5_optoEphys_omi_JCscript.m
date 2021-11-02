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
load Tfig5_imp.mat
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

%% PLOT Average Firing rate from subT n=142 out of 5 sessions 3 mice with at least 7 opto trials  
load('Tfig2_cor.mat');
load SMA_oNO_GoCue142.mat
zSMA_oNO = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);  
load SMA_ocC_GoCue142.mat
zSMA_oCC = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);  

load SMA_omi_GoCue360.mat
% subT_bool360 = logical(sum(Tfig5_imp.nSess == [4,7,10,13,15],2));
SMA = SMA(subT_bool360,:); SSemA = SSemA(subT_bool360,:);
zSMA_omi = (SMA-FRepoch.BLE.mean(subT_bool360))./mean(SSemA,2);
load SMA_cor_GoCue545
SMA = SMA(subT_bool,:); SSemA = SSemA(subT_bool,:);
zSMA_cor = (SMA-FRepoch.BLE.mean(subT_bool))./mean(SSemA,2);

Xt=[1:max(size(nanmean(SMA)))]-parfig.pre;
XoptoStim = Xt(parfig.pre-500: parfig.pre+500);
K=2;
T142=Tfig2_cor(subT_bool,:);
idx_NOL =logical(T142.VMVL & T142.z_exct & ~T142.Opto_inib & ~T142.Opto_exct & T142.Opto_post_sess);
idx_SvTh =logical(T142.VMVL & T142.Opto_inib & T142.z_exct);
idx_SvTh_pos =logical(T142.VMVL & T142.Opto_exct & T142.z_exct);
idx_all = logical(T142.VMVL & T142.z_exct & T142.Opto_post_sess);


%% CORRECT TRIALS vs OPTO COR
%cells types: SvTh- 
%trials: cor vs opto_cor     
close all, 
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
saveas(gcf, 'Average_SvTh-_opto_cor','emf')
%cells types: NOL 
%trials: cor vs opto_cor     
figure,
SMA1=zSMA_cor(idx_NOL,:); col1='k'; idx1=idx_NOL; 
SMA2=zSMA_oCC(idx_NOL,:); col2='b'; idx2=idx_NOL; 
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
saveas(gcf, 'Average_NOL_opto_cor','emf')

%cells types: SvTh+ 
%trials: cor vs opto_cor     
figure,
SMA1=zSMA_cor(idx_SvTh_pos,:); col1='k'; idx1=idx_SvTh_pos; 
SMA2=zSMA_oCC(idx_SvTh_pos,:); col2='m'; idx2=idx_SvTh_pos; 
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
saveas(gcf, 'Average_SvTh+_opto_cor','emf')

%cells types: ALL 
%trials: cor vs opto_cor    
figure,
SMA1=zSMA_cor(idx_all,:); col1='k'; idx1=idx_all; 
SMA2=zSMA_oCC(idx_all,:); col2='y'; idx2=idx_all; 
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
legend(['NONopto cor trials (n=' num2str(sum(idx1)) 'cells ALL)'],[' OPTO cor trials (n=' num2str(sum(idx1)) ' cells ALL)'], 'Location','best')
saveas(gcf, 'Average_ALL_opto_cor','emf')

%% OMI TRIALS vs OPTO OMI
%cells types: SvTh- 
%trials: omi vs opto_omi  
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
saveas(gcf, 'Average_SvTh-_opto_omi','emf')
%
%cells types: NOL 
%trials: omi vs opto_omi  
figure,
SMA1=zSMA_omi(idx_NOL,:); col1='g'; idx1=idx_NOL; 
SMA2=zSMA_oNO(idx_NOL,:); col2='b'; idx2=idx_NOL; 
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
legend(['NONopto omi trials (n=' num2str(sum(idx1)) 'cells NOL)'],[' OPTO omi trials (n=' num2str(sum(idx1)) ' cells NOL)'], 'Location','best')
saveas(gcf, 'Average_NOL_opto_omi','emf')
%cells types: SvTh+ 
%trials: omi vs opto_omi      
figure,
SMA1=zSMA_omi(idx_SvTh_pos,:); col1='g'; idx1=idx_SvTh_pos; 
SMA2=zSMA_oNO(idx_SvTh_pos,:); col2='m'; idx2=idx_SvTh_pos; 
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
legend(['NONopto omi trials (n=' num2str(sum(idx1)) 'cells SvTh+)'],[' OPTO omi trials (n=' num2str(sum(idx1)) ' cells SvTh+)'], 'Location','best')
saveas(gcf, 'Average_SvTh+_opto_omi','emf')
%cells types: ALL 
%trials: omi vs opto_omi    
figure,
SMA1=zSMA_omi(idx_all,:); col1='g'; idx1=idx_all; 
SMA2=zSMA_oNO(idx_all,:); col2='y'; idx2=idx_all; 
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
legend(['NONopto omi trials (n=' num2str(sum(idx1)) 'cells ALL)'],[' OPTO omi trials (n=' num2str(sum(idx1)) ' cells ALL)'], 'Location','best')
saveas(gcf, 'Average_ALL_opto_omi','emf')

%% IMP TRIALS vs OPTO IMP
load SMA_oNO_GoCue142.mat
zSMA_oNO = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);  
load SMA_ocC_GoCue142.mat
zSMA_oCC = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);  
load SMA_imp_GoCue360.mat
% subT_bool360 = logical(sum(Tfig5_imp.nSess == [4,7,10,13,15],2));
SMA = SMA(subT_bool360,:); SSemA = SSemA(subT_bool360,:);
zSMA_imp = (SMA-FRepoch.BLE.mean(subT_bool360))./mean(SSemA,2);
load SMA_cor_GoCue545
SMA = SMA(subT_bool,:); SSemA = SSemA(subT_bool,:);
zSMA_cor = (SMA-FRepoch.BLE.mean(subT_bool))./mean(SSemA,2);

Xt=[1:max(size(nanmean(SMA)))]-parfig.pre;
XoptoStim = Xt(parfig.pre-500: parfig.pre+500);
K=2;
T142=Tfig2_cor(subT_bool,:);
idx_NOL =logical(T142.VMVL & T142.z_exct & ~T142.Opto_inib & ~T142.Opto_exct & T142.Opto_post_sess);
idx_SvTh =logical(T142.VMVL & T142.Opto_inib & T142.z_exct);
idx_SvTh_pos =logical(T142.VMVL & T142.Opto_exct & T142.z_exct);
idx_all = logical(T142.VMVL & T142.z_exct & T142.Opto_post_sess);


%cells types: SvTh- 
%trials: omi vs opto_omi  
figure,
SMA1=zSMA_imp(idx_SvTh,:); col1='g'; idx1=idx_SvTh; 
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
legend(['NONopto imp trials (n=' num2str(sum(idx1)) 'cells SvTh-)'],[' OPTO imp trials (n=' num2str(sum(idx1)) ' cells SvTh-)'], 'Location','best')
saveas(gcf, 'Average_SvTh-_opto_imp','emf')
%
%cells types: NOL 
%trials: imp vs opto_imp  
figure,
SMA1=zSMA_imp(idx_NOL,:); col1='g'; idx1=idx_NOL; 
SMA2=zSMA_oNO(idx_NOL,:); col2='b'; idx2=idx_NOL; 
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
legend(['NONopto imp trials (n=' num2str(sum(idx1)) 'cells NOL)'],[' OPTO imp trials (n=' num2str(sum(idx1)) ' cells NOL)'], 'Location','best')
saveas(gcf, 'Average_NOL_opto_imp','emf')
%cells types: SvTh+ 
%trials: imp vs opto_imp      
figure,
SMA1=zSMA_imp(idx_SvTh_pos,:); col1='g'; idx1=idx_SvTh_pos; 
SMA2=zSMA_oNO(idx_SvTh_pos,:); col2='m'; idx2=idx_SvTh_pos; 
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
legend(['NONopto imp trials (n=' num2str(sum(idx1)) 'cells SvTh+)'],[' OPTO imp trials (n=' num2str(sum(idx1)) ' cells SvTh+)'], 'Location','best')
saveas(gcf, 'Average_SvTh+_opto_imp','emf')
%cells types: ALL 
%trials: imp vs opto_imp    
figure,
SMA1=zSMA_imp(idx_all,:); col1='g'; idx1=idx_all; 
SMA2=zSMA_oNO(idx_all,:); col2='y'; idx2=idx_all; 
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
legend(['NONopto imp trials (n=' num2str(sum(idx1)) 'cells ALL)'],[' OPTO imp trials (n=' num2str(sum(idx1)) ' cells ALL)'], 'Location','best')
saveas(gcf, 'Average_ALL_opto_imp','emf')

