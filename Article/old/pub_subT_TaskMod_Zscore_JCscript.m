% pub_subT_Zscore_CellTypes_TaskUpMod_JCscript
% written by Julien Catanese 03/22/2019

% Set your Threshold for significativity (e.g. 3-std):   
Zthr = 3; 

% load SMA to compute zSMA 
load('SMA_cor_GoCue.mat'); 
GO = parfigUsed.GO; 

% Compute Z-score (zSMA) 
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

zpuf_exct = zSMA_MAX_puf >= Zthr %& Sig_lenght_puf > 75; %at least 50ms above thr 
zdel_exct = zSMA_MAX_del >= Zthr %& Sig_lenght_del > 75;
zres_exct = zSMA_MAX_res >= Zthr %& Sig_lenght_res > 75;

zpuf_inib = zSMA_MIN_puf <= -Zthr ;
zdel_inib = zSMA_MIN_del <= -Zthr ;
zres_inib = zSMA_MIN_res <= -Zthr ;

% VAR to ADD to TABLE
load ('listcell.mat'); T1=[]; T2= []; T1=Tephys; Tephys_z = []; 
T2 = addvars(T1, zSMA_AbsPeak_all, zSMA_MAX_puf, zSMA_MAX_del, zSMA_MAX_res, zSMA_MIN_puf, zSMA_MIN_del, zSMA_MIN_res, zpuf_exct, zdel_exct, zres_exct, zpuf_inib, zdel_inib, zres_inib); 
Tephys_z = T2;  
save('D:\JC_Analysis\listcell.mat','Tephys_z', '-append'); 
Tephys_z(150:155,:) 
disp('Tephys_z SAVED');

%% Count Nb of Cell in each Group: Excit vs Inib vs Complex vs NotSig 
all_excit = Tephys_z.zpuf_exct | Tephys_z.zdel_exct | Tephys_z.zres_exct; 
all_inib  = Tephys_z.zpuf_inib | Tephys_z.zdel_inib | Tephys_z.zres_inib; 

Nexcit_only = sum(Tcombo.VMVL  & all_excit & ~all_inib )
Ninib_only = sum(Tcombo.VMVL  & ~all_excit & all_inib )
Ncomplex_only = sum(Tcombo.VMVL  & all_excit & all_inib )
Nnothing_only = sum(Tcombo.VMVL  & ~all_excit & ~all_inib)

data= [Nexcit_only, Ncomplex_only, Ninib_only]
resolution= 0.05
colmap='bone'
vennX(data, resolution, colmap)
colormapeditor
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

%% SAVE AS A TABLE_Z (TclvR_z)
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

%% Count Nb of Cell in each Group: Excit vs Inib vs Complex vs NotSig 
all_excit=[]; all_inib =[]; 
all_excit = TcLvR_z.cont_zpuf_exct | TcLvR_z.cont_zdel_exct | TcLvR_z.cont_zres_exct; 
all_inib  = TcLvR_z.cont_zpuf_inib | TcLvR_z.cont_zdel_inib | TcLvR_z.cont_zres_inib; 

Nexcit_only_contra = sum(Tcombo.VMVL  & all_excit & ~all_inib )
Ninib_only_contra = sum(Tcombo.VMVL  & ~all_excit & all_inib )
Ncomplex_only_contra = sum(Tcombo.VMVL  & all_excit & all_inib )
Nnothing_only_contra = sum(Tcombo.VMVL  & ~all_excit & ~all_inib)


all_excit=[]; all_inib =[]; 
all_excit = TcLvR_z.ipsi_zpuf_exct | TcLvR_z.ipsi_zdel_exct | TcLvR_z.ipsi_zres_exct; 
all_inib  = TcLvR_z.ipsi_zpuf_inib | TcLvR_z.ipsi_zdel_inib | TcLvR_z.ipsi_zres_inib; 

Nexcit_only_ipsi = sum(Tcombo.VMVL  & all_excit & ~all_inib )
Ninib_only_ipsi = sum(Tcombo.VMVL  & ~all_excit & all_inib )
Ncomplex_only_ipsi = sum(Tcombo.VMVL  & all_excit & all_inib )
Nnothing_only_ipsi = sum(Tcombo.VMVL  & ~all_excit & ~all_inib)

