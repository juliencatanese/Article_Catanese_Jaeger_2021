function [T2]= Update_list_JCfun(T1, SMA, zSMA, FRepoch, ttestEpoch, parfig) 
% function [T2]= Update_list_JCfun(T1, SMA, zSMA, FRepoch, ttestEpoch, parfig) 
% by JC 1/17/2019
T2=[]; 
pre= parfig.pre;

% Epochs mean firing rate: 
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

T2 = addvars(T1(:,1:9), Fr_BLE_Hz, SMA_MAX_Hz, SMA_MIN_Hz, zSMA_AbsPeak_all, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
    zSMA_MAX_puf, zSMA_MAX_del, zSMA_MAX_res, ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P, ... 
    zSMA_MIN_puf, zSMA_MIN_del, zSMA_MIN_res,  ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P); 

row=50
T2(row:row+5,17)
disp('done and done')
Tephys=T2; 
save('D:\JC_Analysis\listcell.mat','Tephys', '-append')

end

