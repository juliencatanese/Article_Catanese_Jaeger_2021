% pub_table3_ephysSUA_task_stat_JCscript
% create and save a table ('Tephys') that contain new cells fields
% Fields:  Fr_BLE_Hz, SMA_CL_MAX_Hz, SMA_CL_MIN_Hz, zSMA_CL_AbsPeak_all, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
%     zSMA_CL_MAX_puf, zSMA_CL_MAX_del, zSMA_CL_MAX_res, ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P, ... 
%     zSMA_CL_MIN_puf, zSMA_CL_MIN_del, zSMA_CL_MIN_res); 
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Tephys')
% written by Julien Catanese 11/26/2018
% last updated JC 2/7/2019.
% last updated JC 3/22/2019: removed the Zscore part. 

%% Create table5 Ipsi-Contra TcLvR (table with  Firing rate info) 
% Ipsi trials (correct Choice Left = cCL)  
ipsi_BLE_Hz = FRepoch_CL.BLE.mean;
ipsi_puf_Hz = FRepoch_CL.puf.mean;
ipsi_del_Hz = FRepoch_CL.del.mean;
ipsi_res_Hz = FRepoch_CL.res.mean;
ipsi_maxSMA_Hz = max(SMA_CL(:,200:end-200)')';
ipsi_minSMA_Hz = min(SMA_CL(:,200:end-200)')';

Hipsi_puf = ttestEpoch_CL.H_puf_BLE;
Hipsi_del = ttestEpoch_CL.H_del_BLE;
Hipsi_res = ttestEpoch_CL.H_res_BLE;

Pipsi_puf = ttestEpoch_CL.P_puf_BLE;
Pipsi_del = ttestEpoch_CL.P_del_BLE;
Pipsi_res = ttestEpoch_CL.P_res_BLE;

% Contra trials (correct Choice Right = cCR)  
cont_BLE_Hz = FRepoch_CR.BLE.mean;
cont_puf_Hz = FRepoch_CR.puf.mean;
cont_del_Hz = FRepoch_CR.del.mean;
cont_res_Hz = FRepoch_CR.res.mean;
cont_maxSMA_Hz = max(SMA_CR(:,200:end-200)')';
cont_minSMA_Hz = min(SMA_CR(:,200:end-200)')';

Hcont_puf = ttestEpoch_CR.H_puf_BLE;
Hcont_del = ttestEpoch_CR.H_del_BLE;
Hcont_res = ttestEpoch_CR.H_res_BLE;

Pcont_puf = ttestEpoch_CR.P_puf_BLE;
Pcont_del = ttestEpoch_CR.P_del_BLE;
Pcont_res = ttestEpoch_CR.P_res_BLE;

% Make table5 = TcLvR
T2=[]; T1=listcell; 
T2 = addvars(T1(:,1:5), ...
cont_BLE_Hz, ipsi_BLE_Hz, Hcont_puf, Hipsi_puf,  Hcont_del, Hipsi_del, Hcont_res, Hipsi_res,...
cont_puf_Hz, ipsi_puf_Hz, cont_del_Hz, ipsi_del_Hz, cont_res_Hz, ipsi_res_Hz,...
cont_maxSMA_Hz, ipsi_maxSMA_Hz, cont_minSMA_Hz, ipsi_minSMA_Hz,... 
Pipsi_puf, Pipsi_del, Pipsi_res, Pcont_puf, Pcont_del, Pcont_res); 
% Save TcLvR into listcell.mat
TcLvR=T2; 
save('D:\JC_Analysis\listcell.mat','TcLvR', '-append'); disp('TcLvR SAVED')
row=250; TcLvR(row:row+15,1:end); disp('TcLvR saved')

%% CONTRA (RIGHT)
% load SMA to compute zSMA 
SMA=[]; zSMA=[]; 
load('SMA_cCR_GoCue.mat'); 
% Compute Z-score (zSMA) 
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2); 
zSMA(find(zSMA==inf))=NaN;
% Get the peack (up or Down) for each epoch of the task
cont_zSMA_MAX_puf = max(zSMA(:,GO-1500:GO-750)')';
cont_zSMA_MAX_del = max(zSMA(:,GO-750:GO)')';
cont_zSMA_MAX_res = max(zSMA(:,GO:GO+750)')';
cont_zSMA_MIN_puf = min(zSMA(:,GO-1500:GO-750)')';
cont_zSMA_MIN_del = min(zSMA(:,GO-750:GO)')';
cont_zSMA_MIN_res = min(zSMA(:,GO:GO+750)')';
cont_zSMA_AbsPeak_all = max(abs(zSMA(:,200:end-200))')';

cont_zpuf_exct = cont_zSMA_MAX_puf >= Zthr ;
cont_zdel_exct = cont_zSMA_MAX_del >= Zthr ;
cont_zres_exct = cont_zSMA_MAX_res >= Zthr ;
cont_zpuf_inib = cont_zSMA_MIN_puf <= -Zthr ;
cont_zdel_inib = cont_zSMA_MIN_del <= -Zthr ;
cont_zres_inib = cont_zSMA_MIN_res <= -Zthr ;
disp('contra done')

%% IPSI (LEFT)
% load SMA to compute zSMA 
SMA=[]; zSMA=[]; 
load('SMA_cCL_GoCue.mat'); 
% Compute Z-score (zSMA) 
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2); 
zSMA(find(zSMA==inf))=NaN;
% Get the peack (up or Down) for each epoch of the task
ipsi_zSMA_MAX_puf = max(zSMA(:,GO-1500:GO-750)')';
ipsi_zSMA_MAX_del = max(zSMA(:,GO-750:GO)')';
ipsi_zSMA_MAX_res = max(zSMA(:,GO:GO+750)')';
ipsi_zSMA_MIN_puf = min(zSMA(:,GO-1500:GO-750)')';
ipsi_zSMA_MIN_del = min(zSMA(:,GO-750:GO)')';
ipsi_zSMA_MIN_res = min(zSMA(:,GO:GO+750)')';
ipsi_zSMA_AbsPeak_all = max(abs(zSMA(:,200:end-200))')';

ipsi_zpuf_exct = ipsi_zSMA_MAX_puf >= Zthr ;
ipsi_zdel_exct = ipsi_zSMA_MAX_del >= Zthr ;
ipsi_zres_exct = ipsi_zSMA_MAX_res >= Zthr ;
ipsi_zpuf_inib = ipsi_zSMA_MIN_puf <= -Zthr ;
ipsi_zdel_inib = ipsi_zSMA_MIN_del <= -Zthr ;
ipsi_zres_inib = ipsi_zSMA_MIN_res <= -Zthr ;
disp('ipsi done')

% SAVE AS A TABLE_Z (TclvR_z)
load ('listcell.mat'); T1=[]; T2= []; T1=TcLvR; TcLvR_z = []; 
T2 = addvars(T1, cont_zSMA_AbsPeak_all, ...
    cont_zSMA_MAX_puf, cont_zSMA_MAX_del, cont_zSMA_MAX_res,...
    cont_zSMA_MIN_puf, cont_zSMA_MIN_del, cont_zSMA_MIN_res,...
    cont_zpuf_exct, cont_zdel_exct, cont_zres_exct, ...
    cont_zpuf_inib, cont_zdel_inib, cont_zres_inib, ...
    ipsi_zSMA_AbsPeak_all, ...
    ipsi_zSMA_MAX_puf, ipsi_zSMA_MAX_del, ipsi_zSMA_MAX_res,...
    ipsi_zSMA_MIN_puf, ipsi_zSMA_MIN_del, ipsi_zSMA_MIN_res,...
    ipsi_zpuf_exct, ipsi_zdel_exct, ipsi_zres_exct, ...
    ipsi_zpuf_inib, ipsi_zdel_inib, ipsi_zres_inib); 

TcLvR_z = T2;  
save('D:\JC_Analysis\listcell.mat','TcLvR_z', '-append'); 
TcLvR_z(150:155,:) 
disp('TcLvR_z SAVED');
