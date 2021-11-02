% pub_table3_ephysSUA_task_stat_JCscript
% create and save a table ('Tephys') that contain new cells fields
% Fields:  Fr_BLE_Hz, SMA_MAX_Hz, SMA_MIN_Hz, zSMA_AbsPeak_all, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
%     zSMA_MAX_puf, zSMA_MAX_del, zSMA_MAX_res, ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P, ... 
%     zSMA_MIN_puf, zSMA_MIN_del, zSMA_MIN_res); 
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Tephys')
% written by Julien Catanese 11/26/2018
% last updated JC 2/7/2019.

%%
clearvars -except parfig nn listcell st av tv


try
    if parfig.saveTABLE ==1
        disp('ATTENTION WILL SAVE A NEW TABLE Tephys')
    end
catch
    parfig.saveTABLE=0
end


% cd('D:\JC_Analysis');
% load('listcell.mat'); 

%% Set parameters for analysis and figures (parfig struct)

% parfig.center_evt =  'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
% parfig.trial_type = {'cor'};
% parfig.col = {'k'};
% parfig.xlim = [-1500 1500];
% parfig.plot=0;

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
% [SMA SSA SSemA FRepoch ttestEpoch trigtimes spxtimes] = listcell_SMA_FRepoch_ttest_JCfun(listcell, parfig);
[SMA FRepoch ttestEpoch trigtimes spxtimes sdf_alltr] = pub_comp_ephySMA_task_ttest_JCfun(listcell, parfig)
% load('SMA_BLE_lastsaved.mat')  %OR  load('SMA_FRepoch_ttest_331cel_cor.mat')

%% zSMA z-score = (x-mean)/(std/sqrt(n))  with x=SMA ; mean=M_BLE and std/sqrt(n) = Sem_BLE
% zSMA = (SMA-FRepoch.BLE.mean)./FRepoch.BLE.Sem;
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2); 
zSMA(find(zSMA==inf))=NaN;

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

zSMA_MAX_puf = max(zSMA(:,GO-1500:GO-750)')';
zSMA_MAX_del = max(zSMA(:,GO-750:GO)')';
zSMA_MAX_res = max(zSMA(:,GO:GO+750)')';
zSMA_MIN_puf = min(zSMA(:,GO-1500:GO-750)')';
zSMA_MIN_del = min(zSMA(:,GO-750:GO)')';
zSMA_MIN_res = min(zSMA(:,GO:GO+750)')';

zSMA_AbsPeak_all = max(abs(zSMA(:,200:end-200))')';

T2 = addvars(T1(:,1:6), SMA_MAX_Hz, SMA_MIN_Hz,  Fr_BLE_Hz, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
     H_puf, H_del, H_res, P_puf, P_del,  P_res, zSMA_MAX_puf, zSMA_MAX_del, zSMA_MAX_res, ... 
    zSMA_MIN_puf, zSMA_MIN_del, zSMA_MIN_res, zSMA_AbsPeak_all); 

Tephys=T2; 
if parfig.saveTABLE ==1 
    save('D:\JC_Analysis\listcell.mat','Tephys', '-append')
    disp('Tephys SAVED')
end

%% display 
row=150;
T2(row:row+5,14:end)
disp('Tephys done')

