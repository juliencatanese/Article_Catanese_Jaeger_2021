% pub_table5_Tfig3_contipsi_dz_JCscript
% create and save a table ('Tephys') that contain new cells fields
% Fields:  Fr_BLE_Hz, SMA_CL_MAX_Hz, SMA_CL_MIN_Hz, zSMA_CL_AbsPeak_all, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
%     zSMA_CL_MAX_puf, zSMA_CL_MAX_del, zSMA_CL_MAX_res, ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P, ... 
%     zSMA_CL_MIN_puf, zSMA_CL_MIN_del, zSMA_CL_MIN_res); 
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Tephys')
% written by Julien Catanese 11/26/2018
% last updated JC 2/7/2019.
% last updated JC 3/22/2019: removed the Zscore part. 
% last updated JC 4/12/2019: readded the Zscore part. save Tfig3_ipsicontra and Tcombo.  

clearvars -except mypath parfig
load('listcell');

%% IPSI: Get Firing Rate Values from FRepoch structure
load('SMA_cCL_GoCue.mat');
SMA_ipsi = SMA; SMA=[]; 
zSMA_ipsi = (SMA_ipsi-FRepoch.BLE.mean)./mean(SSemA,2); 

Fr_BLE_Hz_ipsi = FRepoch.BLE.mean;
Fr_puf_Hz_ipsi = FRepoch.puf.mean;
Fr_del_Hz_ipsi = FRepoch.del.mean;
Fr_res_Hz_ipsi = FRepoch.res.mean;
SMA_MAX_Hz_ipsi = max(SMA_ipsi(:,200:end-200)')';
SMA_MIN_Hz_ipsi = min(SMA_ipsi(:,200:end-200)')';

% IPSI: Get ttest results from ttestEpoch structure
H_puf_ipsi = ttestEpoch.H_puf_BLE;
H_del_ipsi = ttestEpoch.H_del_BLE;
H_res_ipsi = ttestEpoch.H_res_BLE;

P_puf_ipsi = ttestEpoch.P_puf_BLE;
P_del_ipsi = ttestEpoch.P_del_BLE;
P_res_ipsi = ttestEpoch.P_res_BLE;

%% CONTRA: Get Firing Rate Values from FRepoch structure
load('SMA_cCR_GoCue.mat');
SMA_cont = SMA; SMA=[]; 
zSMA_cont = (SMA_cont-FRepoch.BLE.mean)./mean(SSemA,2); 

Fr_BLE_Hz_cont = FRepoch.BLE.mean;
Fr_puf_Hz_cont = FRepoch.puf.mean;
Fr_del_Hz_cont = FRepoch.del.mean;
Fr_res_Hz_cont = FRepoch.res.mean;
SMA_MAX_Hz_cont = max(SMA_cont(:,200:end-200)')';
SMA_MIN_Hz_cont = min(SMA_cont(:,200:end-200)')';

% CONT: Get ttest results from ttestEpoch structure
H_puf_cont = ttestEpoch.H_puf_BLE;
H_del_cont = ttestEpoch.H_del_BLE;
H_res_cont = ttestEpoch.H_res_BLE;

P_puf_cont = ttestEpoch.P_puf_BLE;
P_del_cont = ttestEpoch.P_del_BLE;
P_res_cont = ttestEpoch.P_res_BLE;

%% Z-score ==> Define 7 Classes of cell response during task
Zthr = 3;
dzSMA =zSMA_cont-zSMA_ipsi;

[zbool zvalu] = pub_comp_Zthr_75msEpoch_JCfun(Zthr, dzSMA, parfig)

dzSMA_MAX_puf = zvalu.puf.max; dzSMA_MAX_del = zvalu.del.max; dzSMA_MAX_res = zvalu.res.max;
dzSMA_MIN_puf = zvalu.puf.min; dzSMA_MIN_del = zvalu.del.min; dzSMA_MIN_res = zvalu.res.min;
dzSMA_AbsPeak_all = max(abs(dzSMA(:,200:end-200))')';
dzRL_exct_pufAll = zbool.puf.exct; dzRL_exct_delAll = zbool.del.exct; dzRL_exct_resAll = zbool.res.exct;
dzRL_inib_pufAll = zbool.puf.inib; dzRL_inib_delAll = zbool.del.inib; dzRL_inib_resAll = zbool.res.inib;

dzRL_exctAll      =  dzRL_exct_pufAll | dzRL_exct_delAll | dzRL_exct_resAll;   % All cell that are sig Excited in task
dzRL_inibAll      =  dzRL_inib_pufAll | dzRL_inib_delAll | dzRL_inib_resAll;
dzRL_exct         =  dzRL_exctAll & ~dzRL_inibAll;  % Excite ONLY
dzRL_inib         = ~dzRL_exctAll &  dzRL_inibAll; % Inib ONLY

dzRL_complex      =  dzRL_exctAll &  dzRL_inibAll; % Complex ONLY
dzRL_nosig        = ~dzRL_exctAll & ~dzRL_inibAll; % NoSig ONLY

dzRL_exct_UniMod_1puf      =  dzRL_exct &  dzRL_exct_pufAll & ~dzRL_exct_delAll  & ~dzRL_exct_resAll;
dzRL_exct_UniMod_2del      =  dzRL_exct & ~dzRL_exct_pufAll &  dzRL_exct_delAll  & ~dzRL_exct_resAll;
dzRL_exct_UniMod_3res      =  dzRL_exct & ~dzRL_exct_pufAll & ~dzRL_exct_delAll  &  dzRL_exct_resAll;
dzRL_inib_UniMod_1puf      =  dzRL_inib &  dzRL_inib_pufAll & ~dzRL_inib_delAll  & ~dzRL_inib_resAll;
dzRL_inib_UniMod_2del      =  dzRL_inib & ~dzRL_inib_pufAll &  dzRL_inib_delAll  & ~dzRL_inib_resAll;
dzRL_inib_UniMod_3res      =  dzRL_inib & ~dzRL_inib_pufAll & ~dzRL_inib_delAll  &  dzRL_inib_resAll;

dzRL_exct_biMod_12      =  dzRL_exct &  dzRL_exct_pufAll &  dzRL_exct_delAll  & ~dzRL_exct_resAll;
dzRL_exct_biMod_13      =  dzRL_exct &  dzRL_exct_pufAll & ~dzRL_exct_delAll  &  dzRL_exct_resAll;
dzRL_exct_biMod_23      =  dzRL_exct & ~dzRL_exct_pufAll &  dzRL_exct_delAll  &  dzRL_exct_resAll;
dzRL_inib_biMod_12      =  dzRL_inib &  dzRL_inib_pufAll &  dzRL_inib_delAll  & ~dzRL_inib_resAll;
dzRL_inib_biMod_13      =  dzRL_inib &  dzRL_inib_pufAll & ~dzRL_inib_delAll  &  dzRL_inib_resAll;
dzRL_inib_biMod_23      =  dzRL_inib & ~dzRL_inib_pufAll &  dzRL_inib_delAll  &  dzRL_inib_resAll;

dzRL_exct_triMod_123      =  dzRL_exct &  dzRL_exct_pufAll &  dzRL_exct_delAll  &  dzRL_exct_resAll;
dzRL_inib_triMod_123      =  dzRL_inib &  dzRL_inib_pufAll &  dzRL_inib_delAll  &  dzRL_inib_resAll;

%% SAVING TABLE
if parfig.saveTABLE ==1
    % Tfig3_contipsi_dz = (old TcLvcR)
    Tfig3_contipsi_dz = addvars(Topto, SMA_MAX_Hz_cont, SMA_MIN_Hz_cont, SMA_MAX_Hz_ipsi, SMA_MIN_Hz_ipsi,...
        Fr_BLE_Hz_cont, Fr_puf_Hz_cont, Fr_del_Hz_cont, Fr_res_Hz_cont,...
        Fr_BLE_Hz_ipsi, Fr_puf_Hz_ipsi, Fr_del_Hz_ipsi, Fr_res_Hz_ipsi,...
        dzSMA_AbsPeak_all, dzSMA_MAX_puf, dzSMA_MAX_del,...
        dzSMA_MAX_res, dzSMA_MIN_puf, dzSMA_MIN_del, dzSMA_MIN_res,...
        dzRL_exct_pufAll, dzRL_exct_delAll, dzRL_exct_resAll, ...
        dzRL_inib_pufAll, dzRL_inib_delAll, dzRL_inib_resAll, ...
        dzRL_exctAll, dzRL_inibAll, dzRL_exct, dzRL_inib, dzRL_complex, dzRL_nosig, ...
        dzRL_exct_UniMod_1puf, dzRL_exct_UniMod_2del, dzRL_exct_UniMod_3res, ...
        dzRL_inib_UniMod_1puf, dzRL_inib_UniMod_2del, dzRL_inib_UniMod_3res, ...
        dzRL_exct_biMod_12, dzRL_exct_biMod_13, dzRL_exct_biMod_23, ...
        dzRL_inib_biMod_12, dzRL_inib_biMod_13, dzRL_inib_biMod_23, ...
        dzRL_exct_triMod_123, dzRL_inib_triMod_123,...
        H_puf_cont, H_del_cont, H_res_cont, P_puf_cont, P_del_cont,  P_res_cont,...
        H_puf_ipsi, H_del_ipsi, H_res_ipsi, P_puf_ipsi, P_del_ipsi,  P_res_ipsi);
       
    save('listcell.mat','Tfig3_contipsi_dz', '-append')
    disp('Tfig3_contipsi_dz SAVED'); Tfig3_contipsi_dz(1,:)
        
    % TCOMBO
    Tcombo=[];
    Tcombo=listcell(:,1:5); AreaID=Tcoord.AreaID; VMVL=Tcoord.VMVL;    
    ncell =  Tcoord.ncell; nSess =  Tcoord.nSess; nMouse = Tcoord.nMouse;
    Tcombo = addvars(Tcombo, nMouse, nSess ,VMVL, AreaID);   
    Opto_inib = Topto.Opto_inib; Opto_exct = Topto.Opto_exct; Opto_post_sess= Topto.Opto_post_sess;
    Tcombo= addvars(Tcombo, Opto_inib, Opto_exct, Opto_post_sess);    
    z_exct=Tfig2_cor.z_exct; z_inib=Tfig2_cor.z_inib; z_complex=Tfig2_cor.z_complex;
    z_nosig=Tfig2_cor.z_nosig; z_exct_UniMod_1puf=Tfig2_cor.z_exct_UniMod_1puf;
    z_exct_UniMod_2del=Tfig2_cor.z_exct_UniMod_2del; z_exct_UniMod_3res=Tfig2_cor.z_exct_UniMod_3res;
    z_exct_biMod_12=Tfig2_cor.z_exct_biMod_12; z_exct_biMod_13=Tfig2_cor.z_exct_biMod_13;
    z_exct_biMod_23=Tfig2_cor.z_exct_biMod_23; z_exct_triMod_123=Tfig2_cor.z_exct_triMod_123;   
    Tcombo= addvars(Tcombo, z_exct, z_inib, z_complex, z_nosig, ...
        z_exct_UniMod_1puf, z_exct_UniMod_2del, z_exct_UniMod_3res, ...
        z_exct_biMod_12, z_exct_biMod_13, z_exct_biMod_23, ...
        z_exct_triMod_123);
    dzRL_exct=Tfig3_contipsi_dz.dzRL_exct; 
    dzRL_inib=Tfig3_contipsi_dz.dzRL_inib; 
    dzRL_complex=Tfig3_contipsi_dz.dzRL_complex; 
    dzRL_nosig=Tfig3_contipsi_dz.dzRL_nosig;  
    dzRL_exct_UniMod_1puf=Tfig3_contipsi_dz.dzRL_exct_UniMod_1puf;  
    dzRL_exct_UniMod_2del=Tfig3_contipsi_dz.dzRL_exct_UniMod_2del; 
    dzRL_exct_UniMod_3res=Tfig3_contipsi_dz.dzRL_exct_UniMod_3res; 
    dzRL_exct_biMod_12=Tfig3_contipsi_dz.dzRL_exct_biMod_12; 
    dzRL_exct_biMod_13=Tfig3_contipsi_dz.dzRL_exct_biMod_13; 
    dzRL_exct_biMod_23=Tfig3_contipsi_dz.dzRL_exct_biMod_23; 
    dzRL_exct_triMod_123=Tfig3_contipsi_dz.dzRL_exct_triMod_123;   
    Tcombo= addvars(Tcombo, dzRL_exct, dzRL_inib, dzRL_complex, dzRL_nosig, ...
        dzRL_exct_UniMod_1puf, dzRL_exct_UniMod_2del, dzRL_exct_UniMod_3res, ...
        dzRL_exct_biMod_12, dzRL_exct_biMod_13, dzRL_exct_biMod_23, ...
        dzRL_exct_triMod_123);
    save('Tcombo.mat','Tcombo')
    disp('Tcombo SAVED'); Tcombo(1,:)
    
end


