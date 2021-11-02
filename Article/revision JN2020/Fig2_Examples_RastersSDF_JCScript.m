%  SUA_Article_JCscript
% by JC 11-26-2018
%%
clear all, close all,
mkdir('D:\JC_Figures\Fig_Article');
% cd('D:\JC_Analysis');
cd('D:\DATA EMORY\JC_Analysis')
mypath = 'D:\DATA EMORY\JC_Analysis'; 
% Get_table_listcell_JCcript
load('listcell3.mat');
load('Tfig2_cor.mat')
load('Tcombo.mat')
load('SMA_cor_GoCue545.mat');
parfig.center_evt =  'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.trial_type = {'cor'};
parfig.col = {'k'};
parfig.xlim = [-2750 1250];

% zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2); zSMA(find(zSMA==inf))=NaN;
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSA,2); zSMA(find(zSMA==inf))=NaN;

parfig.WorkFolder = mypath;
parfig.saveSMA=0;
parfig.saveTABLE=0;
parfig.center_evt =  'GoCue';
parfig.BaselineEpoch= [150:2150];
parfig.pre = parfig.BaselineEpoch(end) + 1500 % define how much time before zero
parfig.post = 1500 % define how much time after zero (zero will be defined by trigtimes.cor)
parfig.xlim = [-2750 1250];
parfig.merge=1;

parfig.plotmerge=0;
parfig.nRaster = 2;
% parfig.xlim= [-1750 1000];
parfig.plot=1;
parfig.k=2;
parfig.ylabel = 'z';
parfig.plotshaded= 'sem';

pub_fig2_example_RasterSDF_JCfun(Tfig2_cor, zSMA, listcell, parfig)








%%
%% Figure 2 : Compare Cell firing rate in differents Epochs of the task (Sample vs Delay vs Response)
% Set Parameters for Figure2 = Ephys all correct trials (color Black)
parfig.pre= 2250 + 750 ; pre =parfig.pre ; % 2250ms baseline + 750ms puff
parfig.post= 750 + 900  ; post =parfig.post ;%750ms delay + 750ms response
parfig.BaselineEpoch= [150:2150]; BLE = parfig.BaselineEpoch;
% parfig.BaselineEpoch= [1200:2200]; BLE = parfig.BaselineEpoch;
parfig.center_evt =  'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.trial_type = {'cor'};
parfig.col = {'k'};
parfig.xlim = [-2500 1500]
parfig.plot=0;
SaveON_OFF=0;
%% Fig2A: SUA 3examples (Ephys Epoch)
%% SUA Example#1: Type1 = Sample Epoch
sublist = {'vgat14w14d2S2Ch1clu#01'} %vgat14_w14d8_S1Ch5_clust#1_SDFshaded
cell2plot = listcell(sublist,:);
[SMA1 SSA1 SSemA1 FRepoch1 ttestEpoch1 trigtimes spxtimes sdf_alltr2] = listcell_SMA_FRepoch_ttest_JCfun(cell2plot, parfig)
Plot_SDFshaded_JCfun(SMA1, SSemA1,  cell2plot, parfig, SaveON_OFF)
Plot_Raster_JCfun(SMA1, SSA1, spxtimes, trigtimes, cell2plot, parfig, SaveON_OFF)

%% SUA Example#2: Type2 = Delay Epoch
sublist = {'vgat14w14d2S1Ch6clu#02'}
cell2plot = listcell(sublist,:);
[SMA1 SSA1 SSemA1 FRepoch1 ttestEpoch1 trigtimes spxtimes sdf_alltr] = listcell_SMA_FRepoch_ttest_JCfun(cell2plot, parfig)
Plot_SDFshaded_JCfun(SMA1, SSemA1,  cell2plot, parfig, SaveON_OFF)
% Plot_Raster_JCfun(SMA1, SSA1, spxtimes, trigtimes, cell2plot, parfig, SaveON_OFF)
% Plot_Raster_SDFshaded_JCfun(SMA1, SSemA1, spxtimes, trigtimes, cell2plot, parfig)

%% SUA Example#3 : Type3 = Response Epoch
sublist =  {'vgat17w10d7S4Ch2clu#01'}
cell2plot = listcell(sublist,:);
[SMA1 SSA1 SSemA1 FRepoch1 ttestEpoch1 trigtimes spxtimes sdf_alltr] = listcell_SMA_FRepoch_ttest_JCfun(cell2plot, parfig)
Plot_SDFshaded_JCfun(SMA1, SSemA1,  cell2plot, parfig, SaveON_OFF)

%% SUA Example#5: Type2-3 = Ramp
sublist = {'vgat14w14d8S1Ch5clu#01'}
cell2plot = listcell(sublist,:);
[SMA1 SSA1 SSemA1 FRepoch1 ttestEpoch1 trigtimes spxtimes sdf_alltr2] = listcell_SMA_FRepoch_ttest_JCfun(cell2plot, parfig)
Plot_SDFshaded_JCfun(SMA1, SSemA1,  cell2plot, parfig, SaveON_OFF)
% Plot_Raster_JCfun(SMA1, SSA1, spxtimes, trigtimes, cell2plot, parfig, SaveON_OFF)
% Plot_Raster_SDFshaded_JCfun(SMA1, SSemA1, spxtimes, trigtimes, cell2plot, parfig)
%% SUA Example#6: Type1-2-3 = Ramp
sublist = {'vgat17w10d4S3Ch1clu#01'}
cell2plot = listcell(sublist,:);
[SMA1 SSA1 SSemA1 FRepoch1 ttestEpoch1 trigtimes spxtimes sdf_alltr3] = listcell_SMA_FRepoch_ttest_JCfun(cell2plot, parfig)
Plot_SDFshaded_JCfun(SMA1, SSemA1,  cell2plot, parfig, SaveON_OFF)
% Plot_Raster_JCfun(SMA1, SSA1, spxtimes, trigtimes, cell2plot, parfig, SaveON_OFF)
% Plot_Raster_SDFshaded_JCfun(SMA1, SSemA1, spxtimes, trigtimes, cell2plot, parfig)




%% Fig2 
% Load or Compute SMA and BLE (SdfMeanAll and BaseLineEpoch)
[SMA SSA SSemA FRepoch ttestEpoch trigtimes spxtimes] = listcell_SMA_FRepoch_ttest_JCfun(listcell, parfig)
% load('SMA_FRepoch_ttest_331cel_cor.mat')
%%
load('SMA_BLE_lastsaved.mat')
% plot ALL cell with Delay and Resp ttest P<0.05 ;    
idx_sig_list = find(T2.ttest_del_Hsup & T2.ttest_res_Hsup & T2.zSMA_MAX_del>1)
idx_sig_list = find(T2.ttest_del_Hsup==1 & T2.ttest_res_Hsup==1 & T2.ttest_puf_Hsup==0 & T2.zSMA_MAX_del>3)
close all 
for id=1:size(idx_sig_list,1)
sdf=SMA(idx_sig_list(id),:); shade= SSemA(idx_sig_list(id),:); cell2plot=listcell(idx_sig_list(id),:);
saveON=1;
Plot_SDFshaded_JCfun(sdf, shade, cell2plot, parfig, saveON)
end

