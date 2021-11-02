% pub_table3_ephysSUA_task_stat_JCscript
% create and save a table ('Tephys') that contain new cells fields
% Fields:  Fr_BLE_Hz, SMA_MAX_Hz, SMA_MIN_Hz , Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
%      ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P, ... 
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Tephys')
% written by Julien Catanese 11/26/2018
% last updated JC 2/7/2019.
% last updated JC 3/22/2019: removed the Zscore part. 

%%
clearvars -except mypath parfig nn listcell st av tv

try
    if parfig.saveTABLE ==1
        disp('ATTENTION WILL SAVE A NEW TABLE Tephys')
    end
catch
    parfig.saveTABLE=0
end

%% Set parameters for analysis and figures (parfig struct)
parfig.BaselineEpoch= [150:2150]; 
BLE = parfig.BaselineEpoch;

if parfig.center_evt == 'Delay'
    pre = BLE(end) + 750  ;
    post = 1500 + 750 ;
    GO=pre+750; 
elseif parfig.center_evt == 'GoCue'
    pre = BLE(end) + 1500  ;
    post = 1500 ;
    GO=pre; 
elseif parfig.center_evt == 'Licks'
    pre = BLE(end) + 2250  ;
    post = 1500;
end
parfig.pre=pre;
parfig.post=post;
parfig.GO = GO; 

%% SMA, FrEpoch, ttest (SdfMeanAll and BaseLineEpoch)
[SMA FRepoch ttestEpoch trigtimes spxtimes sdf_alltr] = pub_comp_ephySMA_task_ttest_JCfun(listcell, parfig)

%% Complete list by creating Tephys (table with  Firing rate info) 
T2=[]; T1=listcell; 

Fr_BLE_Hz = FRepoch.BLE.mean;
Fr_puf_Hz = FRepoch.puf.mean;
Fr_del_Hz = FRepoch.del.mean;
Fr_res_Hz = FRepoch.res.mean;
SMA_MAX_Hz = max(SMA(:,200:end-200)')';
SMA_MIN_Hz = min(SMA(:,200:end-200)')';

H_puf = ttestEpoch.H_puf_BLE;
H_del = ttestEpoch.H_del_BLE;
H_res = ttestEpoch.H_res_BLE;

P_puf = ttestEpoch.P_puf_BLE;
P_del = ttestEpoch.P_del_BLE;
P_res = ttestEpoch.P_res_BLE;

T2 = addvars(T1(:,1:6), SMA_MAX_Hz, SMA_MIN_Hz,...
    Fr_BLE_Hz, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
     H_puf, H_del, H_res, P_puf, P_del,  P_res); 

Tephys=T2; 
if parfig.saveTABLE ==1 
    save('D:\JC_Analysis\listcell.mat','Tephys', '-append')
    disp('Tephys SAVED')
end

%% display 
row=150;  T2(row:row+5,14:end); disp('Tephys done')

