% pub_table3_ephysSUA_task_stat_JCscript
% create and save a table ('Tephys') that contain new cells fields
% Fields:  Fr_BLE_Hz, SMA_MAX_Hz, SMA_MIN_Hz, zSMA_AbsPeak_all, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
%     zSMA_MAX_puf, zSMA_MAX_del, zSMA_MAX_res, ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P, ... 
%     zSMA_MIN_puf, zSMA_MIN_del, zSMA_MIN_res); 
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Tephys')
% written by Julien Catanese 11/26/2018
% last updated JC 2/7/2019.

%%
clear all, close all,
mkdir('D:\JC_Figures\Fig_Article');
cd('D:\JC_Analysis');
% Get_table_listcell_JCcript
load('listcell.mat'); 

%% Figure 2 : Compare Cell firing rate in differents Epochs of the task (Sample vs Delay vs Response)
% Set Parameters for Figure2 = Ephys all correct trials (color Black)
parfig.pre= 2250 + 750 ; pre =parfig.pre ; % 2250ms baseline + 750ms puff
parfig.post= 750 + 900  ; post =parfig.post ;%750ms delay + 750ms response
parfig.BaselineEpoch= [150:2150]; BLE = parfig.BaselineEpoch;
% parfig.BaselineEpoch= [1200:2200]; BLE = parfig.BaselineEpoch;
parfig.center_evt =  'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.trial_type = {'cor'};
parfig.col = {'k'};
parfig.xlim = [-1500 1500];
parfig.plot=0;

%% SMA, FrEpoch, ttest (SdfMeanAll and BaseLineEpoch)
[SMA SSA SSemA FRepoch ttestEpoch trigtimes spxtimes] = listcell_SMA_FRepoch_ttest_JCfun(listcell, parfig);
% load('SMA_BLE_lastsaved.mat')  %OR  load('SMA_FRepoch_ttest_331cel_cor.mat')

%% zSMA z-score = (x-mean)/(std/sqrt(n))  with x=SMA ; mean=M_BLE and std/sqrt(n) = Sem_BLE
zSMA = (SMA-FRepoch.BLE.mean)./FRepoch.BLE.Sem;
zSMA(find(zSMA==inf))=NaN;

%% Complete list by creating Tephys (table with  Firing rate info) 
T2=[]; T1=listcell; 

Fr_BLE_Hz = FRepoch.BLE.mean;
Fr_puf_Hz = FRepoch.puf.mean;
Fr_del_Hz = FRepoch.del.mean;
Fr_res_Hz = FRepoch.res.mean;
SMA_MAX_Hz = max(SMA(:,200:end-200)')';
SMA_MIN_Hz = min(SMA(:,200:end-200)')';

ttest_puf_H = ttestEpoch.H_puf_BLE;
ttest_del_H = ttestEpoch.H_del_BLE;
ttest_res_H = ttestEpoch.H_res_BLE;
ttest_puf_P = ttestEpoch.P_puf_BLE;
ttest_del_P = ttestEpoch.P_del_BLE;
ttest_res_P = ttestEpoch.P_res_BLE;


zSMA_MAX_puf = max(zSMA(:,pre-750:pre)')';
zSMA_MAX_del = max(zSMA(:,pre:pre+750)')';
zSMA_MAX_res = max(zSMA(:,pre+750:pre+1500)')';
zSMA_MIN_puf = min(zSMA(:,pre-750:pre)')';
zSMA_MIN_del = min(zSMA(:,pre:pre+750)')';
zSMA_MIN_res = min(zSMA(:,pre+750:pre+1500)')';

zSMA_AbsPeak_all = max(abs(zSMA(:,200:end-200))')';

T2 = addvars(T1(:,1:6), Fr_BLE_Hz, SMA_MAX_Hz, SMA_MIN_Hz, zSMA_AbsPeak_all, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
    zSMA_MAX_puf, zSMA_MAX_del, zSMA_MAX_res, ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P, ... 
    zSMA_MIN_puf, zSMA_MIN_del, zSMA_MIN_res); 

Tephys=T2; 
save('D:\JC_Analysis\listcell.mat','Tephys', '-append')

% display 
row=50;
T2(row:row+5,17)
disp('done and done')

