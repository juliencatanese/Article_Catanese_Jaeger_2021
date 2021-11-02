% pub_table3_ephysSUA_task_stat_JCscript
% create and save a table ('Tephys') that contain new cells fields
% Fields:  Fr_BLE_Hz, SMA_CL_MAX_Hz, SMA_CL_MIN_Hz, zSMA_CL_AbsPeak_all, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
%     zSMA_CL_MAX_puf, zSMA_CL_MAX_del, zSMA_CL_MAX_res, ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P, ... 
%     zSMA_CL_MIN_puf, zSMA_CL_MIN_del, zSMA_CL_MIN_res); 
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Tephys')
% written by Julien Catanese 11/26/2018
% last updated JC 2/7/2019.

%%
clearvars -except parfig nn listcell  st av tv
% cd('D:\JC_Analysis');
load('listcell.mat'); 

%% Set parameters for analysis and figures (parfig struct)
% parfig.center_evt =  'Delay'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
% parfig.trial_type = {'cor'};
% parfig.col = {'k'};
% parfig.xlim = [-1500 1500];
% parfig.plot=0;

parfig.BaselineEpoch= [150:2150]; 
BLE = parfig.BaselineEpoch;
if parfig.center_evt == 'Delay'
    pre = BLE(end) + 750  ;
    post = 1500 + 750 ;
elseif parfig.center_evt == 'GoCue'
    pre = BLE(end) + 1500  ;
    post = 1500 ;
elseif parfig.center_evt == 'Licks'
    pre = BLE(end) + 2250  ;
    post = 1500;
end
parfig.pre=pre;
parfig.post=post;

%% SMA, FrEpoch, ttest (SdfMeanAll and BaseLineEpoch)
% Ipsi trials (correct Choice Left = cCL)  
parfig.trial_type = {'cCL'};
col = {'r'};
[SMA_CL FRepoch_CL ttestEpoch_CL] = pub_comp_ephySMA_task_ttest_JCfun(listcell, parfig);
% Contra trials (correct Choice Right = cCR)  
parfig.trial_type = {'cCR'};
col = {'b'};
[SMA_CR FRepoch_CR ttestEpoch_CR] = pub_comp_ephySMA_task_ttest_JCfun(listcell, parfig);

%% zSMA z-score = (x-mean)/(std/sqrt(n))  with x=SMA ; mean=M_BLE and std/sqrt(n) = Sem_BLE
zSMA_CL = (SMA_CL-FRepoch_CL.BLE.mean)./FRepoch_CL.BLE.Sem;
zSMA_CL(find(zSMA_CL==inf))=NaN;

zSMA_CR = (SMA_CR-FRepoch_CR.BLE.mean)./FRepoch_CR.BLE.Sem;
zSMA_CR(find(zSMA_CR==inf))=NaN;

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

ipsi_maxZSMA_puf_Z = max(zSMA_CL(:,pre-750:pre)')';
ipsi_maxZSMA_del_Z = max(zSMA_CL(:,pre:pre+750)')';
ipsi_maxZSMA_res_Z = max(zSMA_CL(:,pre+750:pre+1500)')';
ipsi_minZSMA_puf_Z = min(zSMA_CL(:,pre-750:pre)')';
ipsi_minZSMA_del_Z = min(zSMA_CL(:,pre:pre+750)')';
ipsi_minZSMA_res_Z = min(zSMA_CL(:,pre+750:pre+1500)')';

ipsi_AbsMaxZSMA_all_Z = max(abs(zSMA_CL(:,200:end-200))')';

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

cont_maxZSMA_puf_Z = max(zSMA_CR(:,pre-750:pre)')';
cont_maxZSMA_del_Z = max(zSMA_CR(:,pre:pre+750)')';
cont_maxZSMA_res_Z = max(zSMA_CR(:,pre+750:pre+1500)')';
cont_minZSMA_puf_Z = min(zSMA_CR(:,pre-750:pre)')';
cont_minZSMA_del_Z = min(zSMA_CR(:,pre:pre+750)')';
cont_minZSMA_res_Z = min(zSMA_CR(:,pre+750:pre+1500)')';

cont_AbsMaxZSMA_all_Z = max(abs(zSMA_CR(:,200:end-200))')';

%% Make table5 = TcLvR
T2=[]; T1=listcell; 
T2 = addvars(T1(:,1:5), ...
cont_BLE_Hz, ipsi_BLE_Hz, Hcont_puf, Hipsi_puf,  Hcont_del, Hipsi_del, Hcont_res, Hipsi_res,...
cont_puf_Hz, ipsi_puf_Hz, cont_del_Hz, ipsi_del_Hz, cont_res_Hz, ipsi_res_Hz,...
cont_maxSMA_Hz, ipsi_maxSMA_Hz, cont_minSMA_Hz, ipsi_minSMA_Hz,...
ipsi_maxZSMA_puf_Z, ipsi_maxZSMA_del_Z, ipsi_maxZSMA_res_Z, ipsi_minZSMA_puf_Z, ipsi_minZSMA_del_Z, ipsi_minZSMA_res_Z,...
cont_maxZSMA_puf_Z, cont_maxZSMA_del_Z, cont_maxZSMA_res_Z, cont_minZSMA_puf_Z, cont_minZSMA_del_Z, cont_minZSMA_res_Z,...
ipsi_AbsMaxZSMA_all_Z, cont_AbsMaxZSMA_all_Z,... 
Pipsi_puf, Pipsi_del, Pipsi_res, Pcont_puf, Pcont_del, Pcont_res); 
%% Save TcLvR into listcell.mat
TcLvR=T2; 
save('D:\JC_Analysis\listcell.mat','TcLvR', '-append');

%% display 
row=250;
TcLvR(row:row+15,1:end);
disp('TcLvR saved')

