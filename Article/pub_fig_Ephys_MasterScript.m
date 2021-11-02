%  SUA_Article_JCscript
% by JC 11-26-2018
%%
clear all, close all,
mkdir('D:\JC_Figures\Fig_Article');
% cd('D:\JC_Analysis');
cd('D:\DATA EMORY\JC_Analysis')
% Get_table_listcell_JCcript
load('listcell3.mat');

%% Figure 2 : Compare Cell firing rate in differents Epochs of the task (Sample vs Delay vs Response)
% Set Parameters for Figure2 = Ephys all correct trials (color Black)
parfig.pre= 2250 + 750 ; pre =parfig.pre ; % 2250ms baseline + 750ms puff
parfig.post= 750 + 900  ; post =parfig.post ;%750ms delay + 750ms response
parfig.BaselineEpoch= [150:2150]; BLE = parfig.BaselineEpoch;
% parfig.BaselineEpoch= [1200:2200]; BLE = parfig.BaselineEpoch;
parfig.center_evt =  'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.trial_type = {'cor'};
parfig.col = {'k'};
parfig.xlim = [-1500 1500]
parfig.plot=0;

%% Fig2A: SUA 3examples (Ephys Epoch)
% SUA Example#1: Type1 = Sample Epoch
% sublist = {'vgat14w14d2z4670S2Ch1clu#1'}
% sublist = {'vgat14w14d8z4900S1Ch5clu#1'} %vgat14_w14d8_S1Ch5_clust#1_SDFshaded
sublist = {'vgat14w14d8S1Ch5clu#01'}
cell2plot = listcell(sublist,:);
[SMA1 SSA1 SSemA1 FRepoch1 ttestEpoch1 trigtimes spxtimes sdf_alltr2] = listcell_SMA_FRepoch_ttest_JCfun(cell2plot, parfig)
%%
SaveON_OFF=0;
Plot_SDFshaded_JCfun(SMA1, SSemA1,  cell2plot, parfig, SaveON_OFF)
Plot_Raster_JCfun(SMA1, SSA1, spxtimes, trigtimes, cell2plot, parfig, SaveON_OFF)
Plot_Raster_SDFshaded_JCfun(SMA1, SSemA1, spxtimes, trigtimes, cell2plot, parfig)
%% SUA Example#6: Type123 = Ramp
sublist = {'vgat17w10d4S3Ch1clu#01'}
cell2plot = listcell(sublist,:);
[SMA1 SSA1 SSemA1 FRepoch1 ttestEpoch1 trigtimes spxtimes sdf_alltr3] = listcell_SMA_FRepoch_ttest_JCfun(cell2plot, parfig)
%%
SaveON_OFF=0;
Plot_SDFshaded_JCfun(SMA1, SSemA1,  cell2plot, parfig, SaveON_OFF)
Plot_Raster_JCfun(SMA1, SSA1, spxtimes, trigtimes, cell2plot, parfig, SaveON_OFF)
% Plot_Raster_SDFshaded_JCfun(SMA1, SSemA1, spxtimes, trigtimes, cell2plot, parfig)
%% SUA Example#2: Type2 = Delay Epoch
sublist = {'vgat14w14d2S1Ch6clu#02'}
cell2plot = listcell(sublist,:);
[SMA1 SSA1 SSemA1 FRepoch1 ttestEpoch1 trigtimes spxtimes sdf_alltr] = listcell_SMA_FRepoch_ttest_JCfun(cell2plot, parfig)
Plot_SDFshaded_JCfun(SMA1, SSemA1,  cell2plot, parfig, SaveON_OFF)
Plot_Raster_JCfun(SMA1, SSA1, spxtimes, trigtimes, cell2plot, parfig, SaveON_OFF)
% Plot_Raster_SDFshaded_JCfun(SMA1, SSemA1, spxtimes, trigtimes, cell2plot, parfig)

%% SUA Example#3 : Type3 = Response Epoch
sublist =  {'vgat17w10d7S4Ch2clu#01'}
cell2plot = listcell(sublist,:);
[SMA1 SSA1 SSemA1 FRepoch1 ttestEpoch1 trigtimes spxtimes sdf_alltr] = listcell_SMA_FRepoch_ttest_JCfun(cell2plot, parfig)
Plot_Raster_SDFshaded_JCfun(SMA1, SSemA1, spxtimes, trigtimes, cell2plot, parfig)

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

%% Fig2B: Population Matrice
% z-score = (x-mean)/(std/sqrt(n))  with x=SMA ; mean=M_BLE and std/sqrt(n) = Sem_BLE
zSMA = (SMA-FRepoch.BLE.mean)./FRepoch.BLE.Sem;
zSMA(find(zSMA==inf))=NaN;

%% Complete the list 
% Call function 
T1=listcell; T2=[]; 
[T2]= Update_list_JCfun(T1, SMA, zSMA, FRepoch, ttestEpoch, parfig);
% A1 = T1(1:5,17) 
% A2= T2(1:5,17)
%% Plot FR and Zscore distrib 
Scatter_Histo_JCfun

%% Select Sig CELLS for pop analysis: 
Zthr= 5
parfig.Zthr=Zthr;

%%%%%%%%%%%%%%%%%%%%
%% 1- Select Based on ZScore
idx_sig_list = find(T2.zSMA_AbsPeak_all > Zthr);
nsc=size(idx_sig_list,1); 
szSMA=zSMA(idx_sig_list,:); 
NszSMA= szSMA./max(abs(szSMA'))';

% Call function to plot Matrices/distrib
%%  (All Epochs)  tt_Sig Normalized Peak disrtrib   
parfig.colormap = 'jet'; 
SortBY = 'peak'; 
Mat2sort=NszSMA;
parfig.title = [SortBY 'Sorted (#Cell='  num2str(nsc) ' Selected by Zhtr >' num2str(Zthr) ' trials-' parfig.trial_type{1} ];
parfig.XsortEpoch =  [pre-1500:pre+post-150];
[sort_Mat2sort] = SortSDFMat_JCfun(Mat2sort, SortBY, 'ascend', parfig);
parfig.caxis = [-1 1]; 
plot_SeqMat_JCfun(sort_Mat2sort, parfig); 
xlim([parfig.XsortEpoch(1) parfig.XsortEpoch(end)]);

%% (Zoom on the GoCue) tt_Sig Normalized Peak disrtrib 
parfig.colormap = 'jet'; 
SortBY = 'peak'
Mat2sort=NszSMA;
parfig.title = [SortBY 'Sorted (#Cell='  num2str(nsc) ' Selected by Zhtr >' num2str(Zthr) ' trials-' parfig.trial_type{1} ];
parfig.XsortEpoch =  [pre:pre+1500];
[sort_Mat2sort] = SortSDFMat_JCfun(Mat2sort, SortBY, 'ascend', parfig);
parfig.caxis = [-1 1]; 
plot_SeqMat_JCfun(sort_Mat2sort, parfig); 
xlim([parfig.XsortEpoch(1) parfig.XsortEpoch(end)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2- Select Based on ttest
idx_sig_list = find(T2.ttest_puf_Hsup==1 | T2.ttest_del_Hsup==1 | T2.ttest_res_Hsup==1);% & T2.zSMA_MAX_del>3);
nsc=size(idx_sig_list,1); 
szSMA=zSMA(idx_sig_list,:); 
NszSMA= szSMA./max(abs(szSMA'))';

% Call function to plot Matrices/distrib
%%  (All Epochs)  tt_Sig Normalized Peak disrtrib   
parfig.colormap = 'jet'; 
SortBY = 'peak'; 
Mat2sort=NszSMA;
parfig.title = [SortBY 'Sorted (#Cell='  num2str(nsc) 'Selected by ttest (any epoch > BLE) trials-' parfig.trial_type{1} ];
parfig.XsortEpoch =  [pre-1500:pre+post-150];
[sort_Mat2sort] = SortSDFMat_JCfun(Mat2sort, SortBY, 'ascend', parfig);
parfig.caxis = [-1 1]; 
plot_SeqMat_JCfun(sort_Mat2sort, parfig); 
xlim([parfig.XsortEpoch(1) parfig.XsortEpoch(end)]);

%% 3- Select Based on Both ttest and ZScore 
%% RAMP CELLS
idx_sig_list = find(T2.ttest_del_Hsup==1 & T2.ttest_res_Hsup==1 & T2.zSMA_MAX_del>5);
nsc=size(idx_sig_list,1); 
szSMA=zSMA(idx_sig_list,:); 
NszSMA= szSMA./max(abs(szSMA'))';

parfig.colormap = 'jet'; 
SortBY = 'peak'
Mat2sort=NszSMA;
parfig.title = [SortBY 'Sorted (#Cell='  num2str(nsc) 'Selected RAMP CELL trials-' parfig.trial_type{1} ];
parfig.XsortEpoch =  [pre:pre+1500];
[sort_Mat2sort] = SortSDFMat_JCfun(Mat2sort, SortBY, 'ascend', parfig);
parfig.caxis = [-1 1]; 
plot_SeqMat_JCfun(sort_Mat2sort, parfig); 
xlim([parfig.XsortEpoch(1) parfig.XsortEpoch(end)]);

%% PUFF CELLS
idx_sig_list = find(T2.ttest_puf_Hsup==1 & T2.zSMA_MAX_del>5);
nsc=size(idx_sig_list,1); 
szSMA=zSMA(idx_sig_list,:); 
NszSMA= szSMA./max(abs(szSMA'))';

parfig.colormap = 'jet'; 
SortBY = 'peak'
Mat2sort=NszSMA;
parfig.title = [SortBY 'Sorted (#Cell='  num2str(nsc) 'Selected PUFF CELL trials-' parfig.trial_type{1} ];
parfig.XsortEpoch =  [pre-1500:pre+750];
[sort_Mat2sort] = SortSDFMat_JCfun(Mat2sort, SortBY, 'ascend', parfig);
parfig.caxis = [-1 1]; 
plot_SeqMat_JCfun(sort_Mat2sort, parfig); 
xlim([parfig.XsortEpoch(1) parfig.XsortEpoch(end)]);

%% PUFF AND RAMP cell 
idx_sig_list = find(T2.ttest_puf_Hsup==1 &  T2.ttest_del_Hsup==1 & T2.ttest_res_Hsup==1 & T2.zSMA_MAX_del>5);
nsc=size(idx_sig_list,1); 
szSMA=zSMA(idx_sig_list,:); 
NszSMA= szSMA./max(abs(szSMA'))';

parfig.colormap = 'jet'; 
SortBY = 'Zthr'
Mat2sort=szSMA;
parfig.title = [SortBY 'Sorted (#Cell='  num2str(nsc) 'Selected PUFF AND RAMP CELL trials-' parfig.trial_type{1} ];
parfig.XsortEpoch =  [pre-100:pre+1000];
[sort_Mat2sort] = SortSDFMat_JCfun(Mat2sort, SortBY, 'ascend', parfig);
parfig.caxis = [-10 25]; 
plot_SeqMat_JCfun(sort_Mat2sort, parfig); parfig.XsortEpoch =  [pre-1500:pre+1500];
xlim([parfig.XsortEpoch(1) parfig.XsortEpoch(end)]);

%% Fig2C: Distribution
figure,
BarSig= zeros(size(zSMA));
ncell=size(BarSig,1);
idxB=find(zSMA>Zthr);
BarSig(idxB)= 1;
parfig.title = 'Distribution'
plot_distrib_3zscore_JCfun(BarSig, ncell, parfig)


%% Figure 3 : Compare Cell responses in Contra vs Ipsi trials
%% Fig3A: 3 SUA examples (Contra vs Ipsi)
close all;
parfig.pre= 2250 + 750 ; pre =parfig.pre ; % 2250ms baseline + 750ms puff
parfig.post= 750 + 900  ; post =parfig.post ;%750ms delay + 750ms response
parfig.BaselineEpoch= [150:2150]; BLE = parfig.BaselineEpoch;
parfig.center_evt =  'Delay'; % center for SDF
parfig.xlim = [-1500 1500]
parfig.plot=0;

% Example Cell #1
sublist =  {'vgat15w10d3z4250S1Ch4clu#1'}
cell2plot = listcell(sublist,:);
parfig.trial_type = {'cCR'};
[SMA1 SSA1 SSemA1 M_BLE1 S_BLE1 Sem_BLE1 trigtimes1 spxtimes1] = listcell_SMA_SdfMeanAll_JCfun(cell2plot, parfig)
parfig.trial_type = {'cCL'};
[SMA2 SSA2 SSemA2 M_BLE2 S_BLE2 Sem_BLE2 trigtimes2 spxtimes2] = listcell_SMA_SdfMeanAll_JCfun(cell2plot, parfig)
parfig.col = [{'b'} {'r'}];
trlim= min(size(trigtimes1,2), size(trigtimes2,2));
Plot_Raster_SDFshaded_JCfun([SMA1;SMA2], [SSemA1;SSemA2] , [spxtimes1; spxtimes2], [trigtimes1(1:trlim);trigtimes2(1:trlim)], cell2plot, parfig)

% Example Cell #2
sublist =  {'vgat14w14d2z4670S2Ch4clu#1'}
cell2plot = listcell(sublist,:);
parfig.trial_type = {'cCR'};
[SMA1 SSA1 SSemA1 M_BLE1 S_BLE1 Sem_BLE1 trigtimes1 spxtimes1] = listcell_SMA_SdfMeanAll_JCfun(cell2plot, parfig)
parfig.trial_type = {'cCL'};
[SMA2 SSA2 SSemA2 M_BLE2 S_BLE2 Sem_BLE2 trigtimes2 spxtimes2] = listcell_SMA_SdfMeanAll_JCfun(cell2plot, parfig)
parfig.col = [{'b'} {'r'}];
trlim= min(size(trigtimes1,2), size(trigtimes2,2));
Plot_Raster_SDFshaded_JCfun([SMA1;SMA2], [SSemA1;SSemA2] , [spxtimes1; spxtimes2], [trigtimes1(1:trlim);trigtimes2(1:trlim)], cell2plot, parfig)

% Example Cell #3
sublist =  {'vgat17w10d7z4300S4Ch2clu#1'}
cell2plot = listcell(sublist,:);
parfig.trial_type = {'cCR'};
[SMA1 SSA1 SSemA1 M_BLE1 S_BLE1 Sem_BLE1 trigtimes1 spxtimes1] = listcell_SMA_SdfMeanAll_JCfun(cell2plot, parfig)
parfig.trial_type = {'cCL'};
[SMA2 SSA2 SSemA2 M_BLE2 S_BLE2 Sem_BLE2 trigtimes2 spxtimes2] = listcell_SMA_SdfMeanAll_JCfun(cell2plot, parfig)
parfig.col = [{'b'} {'r'}];
trlim= min(size(trigtimes1,2), size(trigtimes2,2));
Plot_Raster_SDFshaded_JCfun([SMA1;SMA2], [SSemA1;SSemA2] , [spxtimes1; spxtimes2], [trigtimes1(1:trlim);trigtimes2(1:trlim)], cell2plot, parfig)

%% Fig3B: raw Pop Sequences Matrice (Ipsi vs Contra)
close all;
parfig.pre= 2250 + 750 ; pre =parfig.pre ; % 2250ms baseline + 750ms puff
parfig.post= 750 + 900  ; post =parfig.post ;%750ms delay + 750ms response
parfig.center_evt =  'Delay'; % center for SDF
parfig.xlim = [-1500 1500]
parfig.plot=0;

parfig.trial_type = {'cCL'};
col = {'r'};
[SMA1 SSA1 SSemA1 M_BLE1 S_BLE1 Sem_BLE1] = listcell_SMA_SdfMeanAll_JCfun(listcell, parfig)
save('SMA_BLE_331cel_Left.mat', 'SMA1', 'SSA1', 'SSemA1', 'M_BLE1', 'S_BLE1','Sem_BLE1')

parfig.trial_type = {'cCR'};
col = {'b'};
[SMA2 SSA2 SSemA2 M_BLE2 S_BLE2 Sem_BLE2] = listcell_SMA_SdfMeanAll_JCfun(listcell, parfig)
save('SMA_BLE_331cel_Right.mat', 'SMA2', 'SSA2', 'SSemA2', 'M_BLE2', 'S_BLE2','Sem_BLE2')

%% zscore= (x-mean)/std
BasemeanAll_b = mean(SMA2(:,160:2250),2);
BaseStdAll_b = mean(SSA2(:,160:2250),2);
zSMA2 = (SMA2-BasemeanAll_b)./BaseStdAll_b;

BasemeanAll_r = mean(SMA1(:,160:2250),2);
BaseStdAll_r = mean(SSA2(:,160:2250),2);
zSMA1 = (SMA1-BasemeanAll_r)./BaseStdAll_r;

xx=5
for i=xx:xx+10
    figure(i), plot(zSMA2(i,:),'b'); hold on; plot(zSMA1(i,:),'r')
end
%%  zscore= (x-mean)/std
DzSMA= zSMA1-zSMA2;
BaseDiffmeanAll = mean(DzSMA(:,160:2250),2);
BaseDiffstdAll = std(DzSMA(:,160:2250)')';

zDzSMA = (DzSMA-BaseDiffmeanAll)./BaseDiffstdAll;

xx=5
for i=xx:xx+10
    figure(i), plot(zDzSMA(i,:),'k')%; hold on; plot(zSMA1(i,:),'r')
end
%% Plot Mat for all

NzDzSMA= zDzSMA./max(abs(zDzSMA'))';
Mat2sort=NzDzSMA;
[peak_sorted_mat] = SortPeakSDF_JCfun(Mat2sort, 'ascend', parfig);

parfig.colormap = 'jet'
parfig.caxis = [-1.2 1.2]
parfig.xlim = [pre-1500 pre+1500]
% parfig.title = ['SeqMat-Nsess' num2str(max(size(listcell)) '-#cell' num2str(size(peak_sorted_mat,1)) '-' psth_trig_evt  '-' psth_trial_type{1} ];
plot_SeqMat_JCfun(peak_sorted_mat, parfig);
xlim([pre-1500 pre+post]);

%% Plot Mat for the >3std
Mat2sort =[]; nsig=0;
for i=1:size(zDzSMA,1)
    if ~isempty(find(abs(zDzSMA(i,:))>3))
        Mat2sort = [Mat2sort ; NzDzSMA(i,:)] ;
    else
        nsig=nsig+1;
    end
end
[peak_sorted_mat] = SortPeakSDF_JCfun(Mat2sort, 'ascend', parfig);

parfig.colormap = 'jet'
parfig.caxis = [-1.2 1.2]
parfig.xlim = [pre-1500 pre+1500]
parfig.title = ['SeqMat #cell' num2str(nsig)];
plot_SeqMat_JCfun(peak_sorted_mat, parfig);
xlim([pre-1500 pre+post]);


%% FIGURE 3C DISTRIBUUTION
% plot Dsitribution of Peak increase above 3sem
figure(10),
col = {'r'};
plot_distrib_Peak3sem_JCfun_v2(BarSig, ncell, col, pre, post)

%% Fig3C: Distribution
close all;
figure(10),
BarSig= zeros(size(zDzSMA));
idxB=find(zDzSMA>3);
BarSig(idxB)= 1;
parfig.title = 'Distribution'
parfig.col={'r'}
plot_distrib_3zscore_JCfun(BarSig, nsig, parfig)
parfig.col={'b'}
plot_distrib_3zscore_JCfun(BarSig, nsig, parfig)

%% Fig2C: Distribution
close all; figure(10),
BarSig= zeros(size(zSMA1));
idxB=find(zSMA1>3);
BarSig(idxB)= 1;
parfig.title = 'Distribution'
parfig.col={'b'}
ncell= sum(sum(BarSig,2)>0)
plot_distrib_3zscore_JCfun(BarSig, ncell, parfig)
%%
BarSig= zeros(size(zSMA2));
idxB=find(zSMA2>3);
BarSig(idxB)= 1;
parfig.title = 'Distribution'
parfig.col={'r'}
ncell= sum(sum(BarSig,2)>0)
plot_distrib_3zscore_JCfun(BarSig, ncell, parfig)




%% Figure 4

%% Figure 4 : Compare Cell responses in Impulse vs cor
%% Fig3A: 4 SUA examples (Contra vs Ipsi)
parfig.pre= 2250 + 750 ; pre =parfig.pre ; % 2250ms baseline + 750ms puff
parfig.post= 750 + 900  ; post =parfig.post ;%750ms delay + 750ms response
parfig.center_evt =  'Delay'; % center for SDF
parfig.trial_type = [{'cor'} {'imp'}] % all correct trials
parfig.col=[{'k'},{'g'}]
parfig.xlim = [-1500 1500]
parfig.plot=1
K=1
sublist = {'vgat12w11d5z4300S3Ch2clu#4'}
cell2plot = listcell(sublist,:);
[SMA1 SSA1 SSemA1 M_BLE S_BLE Sem_BLE trigtimes spxtimes] = listcell_SMA_SdfMeanAll_JCfun(cell2plot, parfig)
Plot_Raster_SDFshaded_JCfun(SMA1, SSemA1, spxtimes, trigtimes, cell2plot, parfig)
% SUA Example#2: Type2 = Delay Epoch
sublist = {'vgat12w11d5z4300S3Ch2clu#2'}
cell2plot = listcell(sublist,:);
[SMA1 SSA1 SSemA1 M_BLE S_BLE Sem_BLE trigtimes spxtimes] = listcell_SMA_SdfMeanAll_JCfun(cell2plot, parfig)
Plot_Raster_SDFshaded_JCfun(SMA1, SSemA1, spxtimes, trigtimes, cell2plot, parfig)
% SUA Example#3 : Type3 = Response Epoch
sublist = {'vgat12w11d5z4300S3Ch2clu#3'}
cell2plot = listcell(sublist,:);
[SMA1 SSA1 SSemA1 M_BLE S_BLE Sem_BLE trigtimes spxtimes] = listcell_SMA_SdfMeanAll_JCfun(cell2plot, parfig)
Plot_Raster_SDFshaded_JCfun(SMA1, SSemA1, spxtimes, trigtimes, cell2plot, parfig)

%% Fig4B: raw Pop Sequences Matrice (Ipsi vs Contra)
parfig.pre= 2250 + 750 ; pre =parfig.pre ; % 2250ms baseline + 750ms puff
parfig.post= 750 + 900  ; post =parfig.post ;%750ms delay + 750ms response
parfig.center_evt =  'Delay'; % center for SDF
% parfig.trial_type = [{'cor'} {'imp'}] % all correct trials
% parfig.col=[{'k'},{'g'}]
parfig.xlim = [-1500 1500]
parfig.plot=0
parfig.trial_type = {'cor'};
parfig.col = {'k'};
[SMA1 SSA1 SSemA1] = listcell_SMA_SdfMeanAll_JCfun(listcell, parfig)
parfig.trial_type = {'eNo'};
parfig.col= {'g'};
[SMA2 SSA2 SSemA2] = listcell_SMA_SdfMeanAll_JCfun(listcell, parfig)
%%
close all
% zscore= (x-mean)/std
BasemeanAll_b = mean(SMA2(:,160:2250),2);
BaseStdAll_b = mean(SSA2(:,160:2250),2);
zSMA2 = (SMA2-BasemeanAll_b)./BaseStdAll_b;

BasemeanAll_r = mean(SMA1(:,160:2250),2);
BaseStdAll_r = mean(SSA2(:,160:2250),2);
zSMA1 = (SMA1-BasemeanAll_r)./BaseStdAll_r;

xx=5
for i=xx:xx+10
    figure(i), plot(zSMA2(i,:),'g'); hold on; plot(zSMA1(i,:),'k')
end
%%  zscore= (x-mean)/std
DzSMA= zSMA2-zSMA1;
BaseDiffmeanAll = mean(DzSMA(:,160:2250),2);
BaseDiffstdAll = std(DzSMA(:,160:2250)')';

zDzSMA = (DzSMA-BaseDiffmeanAll)./BaseDiffstdAll;

xx=5
for i=xx:xx+10
    figure(i), plot(zDzSMA(i,:),'r')%; hold on; plot(zSMA1(i,:),'r')
end
%% Plot Mat for all

NzDzSMA= zDzSMA./max(abs(zDzSMA'))';
Mat2sort=NzDzSMA;
[peak_sorted_mat] = SortPeakSDF_JCfun(Mat2sort, 'ascend', parfig);

parfig.colormap = 'cool'
parfig.caxis = [-1 0]
parfig.xlim = [pre-1500 pre+1500]
parfig.title = ['SeqMat-Nsess' 'VGAT12ONLY TEST  -#cell' num2str(size(peak_sorted_mat,1))];
% plot_SeqMat_JCfun(peak_sorted_mat, setplot, pre);
plot_SeqMat_JCfun(peak_sorted_mat, parfig); xlim([pre-1500 pre+post]);


%% Fig4C: Distribution
close all; figure(10),
BarSig= zeros(size(zSMA));
idxB=find(zSMA>K);
BarSig(idxB)= 1;
parfig.title = 'Distribution'
plot_distrib_3zscore_JCfun(BarSig, ncell, parfig)
