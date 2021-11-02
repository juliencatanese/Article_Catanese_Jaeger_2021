% pub_table3_ephysSUA_task_stat_JCscript
% create and save a table ('Tephys') that contain new cells fields
% Fields:  Fr_BLE_Hz, SMA_MAX_Hz, SMA_MIN_Hz , Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
%      ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P, ... 
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Tephys')
% written by Julien Catanese 11/26/2018
% last updated JC 2/7/2019.
% last updated JC 3/22/2019: removed the Zscore part. 

load('SMA_cor_GoCue.mat'); 
%% Create Tephys General with Values 
Tephys=[]; T1=listcell; 

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

Tephys = addvars(T1(:,1:6), SMA_MAX_Hz, SMA_MIN_Hz,...
    Fr_BLE_Hz, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
     H_puf, H_del, H_res, P_puf, P_del,  P_res); 

if parfig.saveTABLE ==1 
    save('D:\JC_Analysis\listcell.mat','Tephys', '-append')
    disp('Tephys SAVED')
end
row=150;  Tephys(row:row+5,14:end); disp('Tephys done')

%% Create Tephys_tt with 7 class of response based on ttest results  
puf_Exc = (Tcoord.VM | Tcoord.VL) &...
    (Tephys.H_puf & (Tephys.Fr_puf_Hz > Tephys.Fr_BLE_Hz))    &...
    (~(Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz))) ; 

del_Exc = (Tcoord.VM | Tcoord.VL) &...
    (Tephys.H_del & (Tephys.Fr_del_Hz > Tephys.Fr_BLE_Hz))    &...
    (~(Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz))) ; 

res_Exc  = (Tcoord.VM | Tcoord.VL) &...
    (Tephys.H_res & (Tephys.Fr_res_Hz > Tephys.Fr_BLE_Hz))    &...
    (~(Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz))) ; 

type7_Exc = (Tcoord.VM | Tcoord.VL) & puf_Exc & del_Exc & res_Exc;
type6_Exc = (Tcoord.VM | Tcoord.VL) & puf_Exc & res_Exc & ~type7_Exc;
type5_Exc = (Tcoord.VM | Tcoord.VL) & del_Exc & res_Exc & ~type7_Exc;
type4_Exc = (Tcoord.VM | Tcoord.VL) & puf_Exc & del_Exc & ~type7_Exc;
type3_Exc = (Tcoord.VM | Tcoord.VL) & res_Exc & ~del_Exc & ~puf_Exc; 
type2_Exc = (Tcoord.VM | Tcoord.VL) & del_Exc & ~res_Exc & ~puf_Exc; 
type1_Exc = (Tcoord.VM | Tcoord.VL) & puf_Exc & ~del_Exc & ~res_Exc; 
type0     = (Tcoord.VM | Tcoord.VL) & Tephys.H_puf==0 &  Tephys.H_del==0 & Tephys.H_res==0; 

Ntt_type0 = sum(type0);
Ntt_type1_Exc = sum(type1_Exc);
Ntt_type2_Exc = sum(type2_Exc);
Ntt_type3_Exc = sum(type3_Exc);
Ntt_type4_Exc = sum(type4_Exc);
Ntt_type5_Exc = sum(type5_Exc);
Ntt_type6_Exc = sum(type6_Exc);
Ntt_type7_Exc = sum(type7_Exc);

Tephys_tt = addvars(T1(:,1:6), SMA_MAX_Hz, SMA_MIN_Hz,...
    Fr_BLE_Hz, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
     H_puf, H_del, H_res, P_puf, P_del,  P_res, ...
     puf_Exc, del_Exc, res_Exc, type0, type1_Exc, type2_Exc, type3_Exc,...
     type4_Exc, type5_Exc, type6_Exc, type7_Exc); 

if parfig.saveTABLE ==1 
    save('D:\JC_Analysis\listcell.mat','Tephys_tt', '-append')
    disp('Tephys SAVED')
end

%% Create Tephys_z with 7 class of response based on zscore results  
Zthr = 3; 
GO = parfig.pre; 
% z-score = (x-mean)/(std/sqrt(n))  with x=SMA ; mean=M_BLE and std/sqrt(n) = SSemA 
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2); zSMA(find(zSMA==inf))=NaN;

% Get the peack (up or Down) for each epoch of the task
zSMA_MAX_puf = max(zSMA(:,GO-1500:GO-750)')';
zSMA_MAX_del = max(zSMA(:,GO-750:GO)')';
zSMA_MAX_res = max(zSMA(:,GO:GO+750)')';
zSMA_MIN_puf = min(zSMA(:,GO-1500:GO-750)')';
zSMA_MIN_del = min(zSMA(:,GO-750:GO)')';
zSMA_MIN_res = min(zSMA(:,GO:GO+750)')';
zSMA_AbsPeak_all = max(abs(zSMA(:,200:end-200))')';

% Sigexct_lenght_puf =  sum(zSMA(:,GO-1500:GO-750)> Zthr, 2)
% Sigexct_lenght_del =  sum(zSMA(:,GO-750:GO)> Zthr, 2)
% Sigexct_lenght_res =  sum(zSMA(:,GO:GO+750)> Zthr, 2)
% Siginib_lenght_puf =  sum(zSMA(:,GO-1500:GO-750)< Zthr, 2)
% Siginib_lenght_del =  sum(zSMA(:,GO-750:GO)< Zthr, 2)
% Siginib_lenght_res =  sum(zSMA(:,GO:GO+750)< Zthr, 2)

zpuf_exct = zSMA_MAX_puf >= Zthr; %& Sig_lenght_puf > 75; %at least 50ms above thr 
zdel_exct = zSMA_MAX_del >= Zthr; %& Sig_lenght_del > 75;
zres_exct = zSMA_MAX_res >= Zthr; %& Sig_lenght_res > 75;

zpuf_inib = zSMA_MIN_puf <= -Zthr ;
zdel_inib = zSMA_MIN_del <= -Zthr ;
zres_inib = zSMA_MIN_res <= -Zthr ;

% VAR to ADD to TABLE
load ('listcell.mat'); T1=[];  T1=Tephys; Tephys_z = []; 
Tephys_z = addvars(T1, zSMA_AbsPeak_all, zSMA_MAX_puf, zSMA_MAX_del,...
    zSMA_MAX_res, zSMA_MIN_puf, zSMA_MIN_del, zSMA_MIN_res, zpuf_exct,...
    zdel_exct, zres_exct, zpuf_inib, zdel_inib, zres_inib); 
save('D:\JC_Analysis\listcell.mat','Tephys_z', '-append'); 
Tephys_z(150:155,:) 
disp('Tephys_z SAVED');


% all_excit = Tephys_z.zpuf_exct | Tephys_z.zdel_exct | Tephys_z.zres_exct; 
% all_inib  = Tephys_z.zpuf_inib | Tephys_z.zdel_inib | Tephys_z.zres_inib; 
% 
% Nexcit_only = sum(Tcombo.VMVL  & all_excit & ~all_inib )
% Ninib_only = sum(Tcombo.VMVL  & ~all_excit & all_inib )
% Ncomplex_only = sum(Tcombo.VMVL  & all_excit & all_inib )
% Nnothing_only = sum(Tcombo.VMVL  & ~all_excit & ~all_inib)
% 
% data= [Nexcit_only, Ncomplex_only, Ninib_only]
% resolution= 0.05
% colmap='bone'
% vennX(data, resolution, colmap)
% colormapeditor