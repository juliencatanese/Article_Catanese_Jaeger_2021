% pub_table6_Tfig5_imp_JCscript
% create and save a table ('Tephys') that omiain new cells fields
% Fields:  Fr_BLE_Hz, SMA_CL_MAX_Hz, SMA_CL_MIN_Hz, zSMA_CL_AbsPeak_all, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
%     zSMA_CL_MAX_puf, zSMA_CL_MAX_del, zSMA_CL_MAX_res, ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P, ...
%     zSMA_CL_MIN_puf, zSMA_CL_MIN_del, zSMA_CL_MIN_res);
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Tephys')
% written by Julien Catanese 11/26/2018
% last updated JC 2/7/2019.
% last updated JC 3/22/2019: removed the Zscore part.
% last updated JC 4/12/2019: readded the Zscore part. save Tfig3_impomira and Tcombo.

% clearvars -except mypath parfig sub_listcell
%% Get zSMA_imp: Get Firing Rate Values from FRepoch structure
load('SMA_imp_GoCue360.mat');
SMA_imp = SMA; SMA=[];
zSMA_imp = (SMA_imp-FRepoch.BLE.mean)./mean(SSemA,2);

Fr_BLE_Hz_imp = FRepoch.BLE.mean; Fr_puf_Hz_imp = FRepoch.puf.mean;
Fr_del_Hz_imp = FRepoch.del.mean; Fr_res_Hz_imp = FRepoch.res.mean;
SMA_MAX_Hz_imp = max(SMA_imp(:,200:end-200)')'; SMA_MIN_Hz_imp = min(SMA_imp(:,200:end-200)')';
H_puf_imp = ttestEpoch.H_puf_BLE; P_puf_imp = ttestEpoch.P_puf_BLE;
H_del_imp = ttestEpoch.H_del_BLE; P_del_imp = ttestEpoch.P_del_BLE;
H_res_imp = ttestEpoch.H_res_BLE; P_res_imp = ttestEpoch.P_res_BLE;

%% Get zSMA_omi: Get Firing Rate Values from FRepoch structure
load('SMA_omi_GoCue360.mat');
SMA_omi = SMA; SMA=[];
zSMA_omi = (SMA_omi-FRepoch.BLE.mean)./mean(SSemA,2);

Fr_BLE_Hz_omi = FRepoch.BLE.mean; Fr_puf_Hz_omi = FRepoch.puf.mean;
Fr_del_Hz_omi = FRepoch.del.mean; Fr_res_Hz_omi = FRepoch.res.mean;
SMA_MAX_Hz_omi = max(SMA_omi(:,200:end-200)')'; SMA_MIN_Hz_omi = min(SMA_omi(:,200:end-200)')';
H_puf_omi = ttestEpoch.H_puf_BLE; P_puf_omi = ttestEpoch.P_puf_BLE;
H_del_omi = ttestEpoch.H_del_BLE; P_del_omi = ttestEpoch.P_del_BLE;
H_res_omi = ttestEpoch.H_res_BLE; P_res_omi = ttestEpoch.P_res_BLE;

%% Get zSMA_cor
load('SMA_cor_GoCue545.mat');
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
zSMA_cor= zSMA(sub_listcell.ncell,:);
%% Z-score ==> Define 7 Classes of cell response during task
dzSMA_imp =zSMA_imp-zSMA_cor;
dzSMA_omi =zSMA_omi-zSMA_cor;

%% IMPULSE TRIALS
Zthr = 3 % was 3  
[zbool zvalu] = pub_comp_Zthr_75msEpoch_JCfun(Zthr, dzSMA_imp, parfig)

dzicSMA_MAX_puf = zvalu.puf.max; dzicSMA_MAX_del = zvalu.del.max; dzicSMA_MAX_res = zvalu.res.max;
dzicSMA_MIN_puf = zvalu.puf.min; dzicSMA_MIN_del = zvalu.del.min; dzicSMA_MIN_res = zvalu.res.min;
dzicSMA_AbsPeak_all = max(abs(dzSMA_imp(:,200:end-200))')';

dzic_exct_pufAll = zbool.puf.exct; dzic_exct_delAll = zbool.del.exct; dzic_exct_resAll = zbool.res.exct;
dzic_inib_pufAll = zbool.puf.inib; dzic_inib_delAll = zbool.del.inib; dzic_inib_resAll = zbool.res.inib;

dzic_exctAll      =  dzic_exct_pufAll | dzic_exct_delAll | dzic_exct_resAll;   % All cell that are sig Excited in task
dzic_inibAll      =  dzic_inib_pufAll | dzic_inib_delAll | dzic_inib_resAll;
dzic_exct         =  dzic_exctAll & ~dzic_inibAll;  % Excite ONLY
dzic_inib         = ~dzic_exctAll &  dzic_inibAll; % Inib ONLY

dzic_complex      =  dzic_exctAll &  dzic_inibAll; % Complex ONLY
dzic_nosig        = ~dzic_exctAll & ~dzic_inibAll; % NoSig ONLY

dzic_exct_UniMod_1puf      =  dzic_exct &  dzic_exct_pufAll & ~dzic_exct_delAll  & ~dzic_exct_resAll;
dzic_exct_UniMod_2del      =  dzic_exct & ~dzic_exct_pufAll &  dzic_exct_delAll  & ~dzic_exct_resAll;
dzic_exct_UniMod_3res      =  dzic_exct & ~dzic_exct_pufAll & ~dzic_exct_delAll  &  dzic_exct_resAll;
dzic_inib_UniMod_1puf      =  dzic_inib &  dzic_inib_pufAll & ~dzic_inib_delAll  & ~dzic_inib_resAll;
dzic_inib_UniMod_2del      =  dzic_inib & ~dzic_inib_pufAll &  dzic_inib_delAll  & ~dzic_inib_resAll;
dzic_inib_UniMod_3res      =  dzic_inib & ~dzic_inib_pufAll & ~dzic_inib_delAll  &  dzic_inib_resAll;

dzic_exct_biMod_12      =  dzic_exct &  dzic_exct_pufAll &  dzic_exct_delAll  & ~dzic_exct_resAll;
dzic_exct_biMod_13      =  dzic_exct &  dzic_exct_pufAll & ~dzic_exct_delAll  &  dzic_exct_resAll;
dzic_exct_biMod_23      =  dzic_exct & ~dzic_exct_pufAll &  dzic_exct_delAll  &  dzic_exct_resAll;
dzic_inib_biMod_12      =  dzic_inib &  dzic_inib_pufAll &  dzic_inib_delAll  & ~dzic_inib_resAll;
dzic_inib_biMod_13      =  dzic_inib &  dzic_inib_pufAll & ~dzic_inib_delAll  &  dzic_inib_resAll;
dzic_inib_biMod_23      =  dzic_inib & ~dzic_inib_pufAll &  dzic_inib_delAll  &  dzic_inib_resAll;

dzic_exct_triMod_123      =  dzic_exct &  dzic_exct_pufAll &  dzic_exct_delAll  &  dzic_exct_resAll;
dzic_inib_triMod_123      =  dzic_inib &  dzic_inib_pufAll &  dzic_inib_delAll  &  dzic_inib_resAll;


%% OMISSION TRIALS
[zbool zvalu] = pub_comp_Zthr_75msEpoch_JCfun(Zthr, dzSMA_omi, parfig)

dzocSMA_MAX_puf = zvalu.puf.max; dzocSMA_MAX_del = zvalu.del.max; dzocSMA_MAX_res = zvalu.res.max;
dzocSMA_MIN_puf = zvalu.puf.min; dzocSMA_MIN_del = zvalu.del.min; dzocSMA_MIN_res = zvalu.res.min;
dzocSMA_AbsPeak_all = max(abs(dzSMA_omi(:,200:end-200))')';
dzoc_exct_pufAll = zbool.puf.exct; dzoc_exct_delAll = zbool.del.exct; dzoc_exct_resAll = zbool.res.exct;
dzoc_inib_pufAll = zbool.puf.inib; dzoc_inib_delAll = zbool.del.inib; dzoc_inib_resAll = zbool.res.inib;

dzoc_exctAll      =  dzoc_exct_pufAll | dzoc_exct_delAll | dzoc_exct_resAll;   % All cell that are sig Excited in task
dzoc_inibAll      =  dzoc_inib_pufAll | dzoc_inib_delAll | dzoc_inib_resAll;
dzoc_exct         =  dzoc_exctAll & ~dzoc_inibAll;  % Excite ONLY
dzoc_inib         = ~dzoc_exctAll &  dzoc_inibAll; % Inib ONLY

dzoc_complex      =  dzoc_exctAll &  dzoc_inibAll; % Complex ONLY
dzoc_nosig        = ~dzoc_exctAll & ~dzoc_inibAll; % NoSig ONLY

dzoc_exct_UniMod_1puf      =  dzoc_exct &  dzoc_exct_pufAll & ~dzoc_exct_delAll  & ~dzoc_exct_resAll;
dzoc_exct_UniMod_2del      =  dzoc_exct & ~dzoc_exct_pufAll &  dzoc_exct_delAll  & ~dzoc_exct_resAll;
dzoc_exct_UniMod_3res      =  dzoc_exct & ~dzoc_exct_pufAll & ~dzoc_exct_delAll  &  dzoc_exct_resAll;
dzoc_inib_UniMod_1puf      =  dzoc_inib &  dzoc_inib_pufAll & ~dzoc_inib_delAll  & ~dzoc_inib_resAll;
dzoc_inib_UniMod_2del      =  dzoc_inib & ~dzoc_inib_pufAll &  dzoc_inib_delAll  & ~dzoc_inib_resAll;
dzoc_inib_UniMod_3res      =  dzoc_inib & ~dzoc_inib_pufAll & ~dzoc_inib_delAll  &  dzoc_inib_resAll;

dzoc_exct_biMod_12      =  dzoc_exct &  dzoc_exct_pufAll &  dzoc_exct_delAll  & ~dzoc_exct_resAll;
dzoc_exct_biMod_13      =  dzoc_exct &  dzoc_exct_pufAll & ~dzoc_exct_delAll  &  dzoc_exct_resAll;
dzoc_exct_biMod_23      =  dzoc_exct & ~dzoc_exct_pufAll &  dzoc_exct_delAll  &  dzoc_exct_resAll;
dzoc_inib_biMod_12      =  dzoc_inib &  dzoc_inib_pufAll &  dzoc_inib_delAll  & ~dzoc_inib_resAll;
dzoc_inib_biMod_13      =  dzoc_inib &  dzoc_inib_pufAll & ~dzoc_inib_delAll  &  dzoc_inib_resAll;
dzoc_inib_biMod_23      =  dzoc_inib & ~dzoc_inib_pufAll &  dzoc_inib_delAll  &  dzoc_inib_resAll;

dzoc_exct_triMod_123      =  dzoc_exct &  dzoc_exct_pufAll &  dzoc_exct_delAll  &  dzoc_exct_resAll;
dzoc_inib_triMod_123      =  dzoc_inib &  dzoc_inib_pufAll &  dzoc_inib_delAll  &  dzoc_inib_resAll;


%% SAVING TABLE
if parfig.saveTABLE ==1
    T1 = []; T1=removevars(Tfig1_VMopto,[1,10,11,12,18,19,20,21,22 ])
    % Tfig5_imp = (old TcLvcR)
    Tfig5_imp = addvars(T1(sub_listcell.ncell,:), SMA_MAX_Hz_omi, SMA_MIN_Hz_omi, SMA_MAX_Hz_imp, SMA_MIN_Hz_imp,...
        Fr_BLE_Hz_omi, Fr_puf_Hz_omi, Fr_del_Hz_omi, Fr_res_Hz_omi,...
        Fr_BLE_Hz_imp, Fr_puf_Hz_imp, Fr_del_Hz_imp, Fr_res_Hz_imp,...
        dzicSMA_AbsPeak_all, dzicSMA_MAX_puf, dzicSMA_MAX_del,...
        dzicSMA_MAX_res, dzicSMA_MIN_puf, dzicSMA_MIN_del, dzicSMA_MIN_res,...
        dzic_exct_pufAll, dzic_exct_delAll, dzic_exct_resAll, ...
        dzic_inib_pufAll, dzic_inib_delAll, dzic_inib_resAll, ...
        dzic_exctAll, dzic_inibAll, dzic_exct, dzic_inib, dzic_complex, dzic_nosig, ...
        dzic_exct_UniMod_1puf, dzic_exct_UniMod_2del, dzic_exct_UniMod_3res, ...
        dzic_inib_UniMod_1puf, dzic_inib_UniMod_2del, dzic_inib_UniMod_3res, ...
        dzic_exct_biMod_12, dzic_exct_biMod_13, dzic_exct_biMod_23, ...
        dzic_inib_biMod_12, dzic_inib_biMod_13, dzic_inib_biMod_23, ...
        dzic_exct_triMod_123, dzic_inib_triMod_123,...
        dzocSMA_AbsPeak_all, dzocSMA_MAX_puf, dzocSMA_MAX_del,...
        dzocSMA_MAX_res, dzocSMA_MIN_puf, dzocSMA_MIN_del, dzocSMA_MIN_res,...
        dzoc_exct_pufAll, dzoc_exct_delAll, dzoc_exct_resAll, ...
        dzoc_inib_pufAll, dzoc_inib_delAll, dzoc_inib_resAll, ...
        dzoc_exctAll, dzoc_inibAll, dzoc_exct, dzoc_inib, dzoc_complex, dzoc_nosig, ...
        dzoc_exct_UniMod_1puf, dzoc_exct_UniMod_2del, dzoc_exct_UniMod_3res, ...
        dzoc_inib_UniMod_1puf, dzoc_inib_UniMod_2del, dzoc_inib_UniMod_3res, ...
        dzoc_exct_biMod_12, dzoc_exct_biMod_13, dzoc_exct_biMod_23, ...
        dzoc_inib_biMod_12, dzoc_inib_biMod_13, dzoc_inib_biMod_23, ...
        dzoc_exct_triMod_123, dzoc_inib_triMod_123,...
        H_puf_omi, H_del_omi, H_res_omi, P_puf_omi, P_del_omi,  P_res_omi,...
        H_puf_imp, H_del_imp, H_res_imp, P_puf_imp, P_del_imp,  P_res_imp);
    
%     save('listcell.mat','Tfig5_imp', '-append')
    disp('Tfig5_imp SAVED'); Tfig5_imp(1,:)
    
    
end


