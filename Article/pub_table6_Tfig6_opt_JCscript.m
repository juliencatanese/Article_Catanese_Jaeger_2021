% pub_table6_Tfig5_opt_JCscript
% create and save a table ('Tephys') that omiain new cells fields
% Fields:  Fr_BLE_Hz, SMA_CL_MAX_Hz, SMA_CL_MIN_Hz, zSMA_CL_AbsPeak_all, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
%     zSMA_CL_MAX_puf, zSMA_CL_MAX_del, zSMA_CL_MAX_res, ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P, ...
%     zSMA_CL_MIN_puf, zSMA_CL_MIN_del, zSMA_CL_MIN_res);
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Tephys')
% written by Julien Catanese 11/26/2018
% last updated JC 2/7/2019.
% last updated JC 3/22/2019: removed the Zscore part.
% last updated JC 4/12/2019: readded the Zscore part. save Tfig3_optomira and Tcombo.

% clearvars -except mypath parfig sub_listcell
%% Get zSMA_opt: Get Firing Rate Values from FRepoch structure
SMA_opt = SMA; SMA=[];
zSMA_opt = (SMA_opt-FRepoch.BLE.mean)./mean(SSemA,2);

Fr_BLE_Hz_opt = FRepoch.BLE.mean; Fr_puf_Hz_opt = FRepoch.puf.mean;
Fr_del_Hz_opt = FRepoch.del.mean; Fr_res_Hz_opt = FRepoch.res.mean;
SMA_MAX_Hz_opt = max(SMA_opt(:,200:end-200)')'; SMA_MIN_Hz_opt = min(SMA_opt(:,200:end-200)')';
H_puf_opt = ttestEpoch.H_puf_BLE; P_puf_opt = ttestEpoch.P_puf_BLE;
H_del_opt = ttestEpoch.H_del_BLE; P_del_opt = ttestEpoch.P_del_BLE;
H_res_opt = ttestEpoch.H_res_BLE; P_res_opt = ttestEpoch.P_res_BLE;

%% Get zSMA_cor
load('SMA_cor_GoCue545.mat');
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
zSMA_cor= zSMA(sub_listcell.ncell,:);

%% Z-score ==> Define 7 Classes of cell response during task
dzSMA_opt =zSMA_opt-zSMA_cor;

%% DZ_opto table 
Zthr = 3;
[zbool zvalu] = pub_comp_Zthr_75msEpoch_JCfun(Zthr, dzSMA_opt, parfig)

dzoptSMA_MAX_puf = zvalu.puf.max; dzoptSMA_MAX_del = zvalu.del.max; dzoptSMA_MAX_res = zvalu.res.max;
dzoptSMA_MIN_puf = zvalu.puf.min; dzoptSMA_MIN_del = zvalu.del.min; dzoptSMA_MIN_res = zvalu.res.min;
dzoptSMA_AbsPeak_all = max(abs(dzSMA_opt(:,200:end-200))')';

dzopt_exct_pufAll = zbool.puf.exct; dzopt_exct_delAll = zbool.del.exct; dzopt_exct_resAll = zbool.res.exct;
dzopt_inib_pufAll = zbool.puf.inib; dzopt_inib_delAll = zbool.del.inib; dzopt_inib_resAll = zbool.res.inib;

dzopt_exctAll      =  dzopt_exct_pufAll | dzopt_exct_delAll | dzopt_exct_resAll;   % All cell that are sig Excited in task
dzopt_inibAll      =  dzopt_inib_pufAll | dzopt_inib_delAll | dzopt_inib_resAll;
dzopt_exct         =  dzopt_exctAll & ~dzopt_inibAll;  % Excite ONLY
dzopt_inib         = ~dzopt_exctAll &  dzopt_inibAll; % Inib ONLY

dzopt_complex      =  dzopt_exctAll &  dzopt_inibAll; % Complex ONLY
dzopt_nosig        = ~dzopt_exctAll & ~dzopt_inibAll; % NoSig ONLY

dzopt_exct_UniMod_1puf      =  dzopt_exct &  dzopt_exct_pufAll & ~dzopt_exct_delAll  & ~dzopt_exct_resAll;
dzopt_exct_UniMod_2del      =  dzopt_exct & ~dzopt_exct_pufAll &  dzopt_exct_delAll  & ~dzopt_exct_resAll;
dzopt_exct_UniMod_3res      =  dzopt_exct & ~dzopt_exct_pufAll & ~dzopt_exct_delAll  &  dzopt_exct_resAll;
dzopt_inib_UniMod_1puf      =  dzopt_inib &  dzopt_inib_pufAll & ~dzopt_inib_delAll  & ~dzopt_inib_resAll;
dzopt_inib_UniMod_2del      =  dzopt_inib & ~dzopt_inib_pufAll &  dzopt_inib_delAll  & ~dzopt_inib_resAll;
dzopt_inib_UniMod_3res      =  dzopt_inib & ~dzopt_inib_pufAll & ~dzopt_inib_delAll  &  dzopt_inib_resAll;

dzopt_exct_biMod_12      =  dzopt_exct &  dzopt_exct_pufAll &  dzopt_exct_delAll  & ~dzopt_exct_resAll;
dzopt_exct_biMod_13      =  dzopt_exct &  dzopt_exct_pufAll & ~dzopt_exct_delAll  &  dzopt_exct_resAll;
dzopt_exct_biMod_23      =  dzopt_exct & ~dzopt_exct_pufAll &  dzopt_exct_delAll  &  dzopt_exct_resAll;
dzopt_inib_biMod_12      =  dzopt_inib &  dzopt_inib_pufAll &  dzopt_inib_delAll  & ~dzopt_inib_resAll;
dzopt_inib_biMod_13      =  dzopt_inib &  dzopt_inib_pufAll & ~dzopt_inib_delAll  &  dzopt_inib_resAll;
dzopt_inib_biMod_23      =  dzopt_inib & ~dzopt_inib_pufAll &  dzopt_inib_delAll  &  dzopt_inib_resAll;

dzopt_exct_triMod_123      =  dzopt_exct &  dzopt_exct_pufAll &  dzopt_exct_delAll  &  dzopt_exct_resAll;
dzopt_inib_triMod_123      =  dzopt_inib &  dzopt_inib_pufAll &  dzopt_inib_delAll  &  dzopt_inib_resAll;


%% SAVING TABLE
if parfig.saveTABLE ==1
    T1 = []; T1=removevars(Tfig1_VMopto,[1,10,11,12,18,19,20,21,22 ])
    % Tfig6_opt = (old TcLvcR)
    Tfig6_opt = addvars(T1(sub_listcell.ncell,:), SMA_MAX_Hz_opt, SMA_MIN_Hz_opt,...
        Fr_BLE_Hz_opt, Fr_puf_Hz_opt, Fr_del_Hz_opt, Fr_res_Hz_opt,...
        dzoptSMA_AbsPeak_all, dzoptSMA_MAX_puf, dzoptSMA_MAX_del,...
        dzoptSMA_MAX_res, dzoptSMA_MIN_puf, dzoptSMA_MIN_del, dzoptSMA_MIN_res,...
        dzopt_exct_pufAll, dzopt_exct_delAll, dzopt_exct_resAll, ...
        dzopt_inib_pufAll, dzopt_inib_delAll, dzopt_inib_resAll, ...
        dzopt_exctAll, dzopt_inibAll, dzopt_exct, dzopt_inib, dzopt_complex, dzopt_nosig, ...
        dzopt_exct_UniMod_1puf, dzopt_exct_UniMod_2del, dzopt_exct_UniMod_3res, ...
        dzopt_inib_UniMod_1puf, dzopt_inib_UniMod_2del, dzopt_inib_UniMod_3res, ...
        dzopt_exct_biMod_12, dzopt_exct_biMod_13, dzopt_exct_biMod_23, ...
        dzopt_inib_biMod_12, dzopt_inib_biMod_13, dzopt_inib_biMod_23, ...
        dzopt_exct_triMod_123, dzopt_inib_triMod_123,...
        H_puf_opt, H_del_opt, H_res_opt, P_puf_opt, P_del_opt,  P_res_opt);
    
    save('listcell.mat','Tfig6_opt', '-append')
    disp('Tfig6_opt SAVED'); Tfig6_opt(1,:)
end


