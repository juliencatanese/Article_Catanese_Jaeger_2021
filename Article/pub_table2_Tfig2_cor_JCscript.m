% pub_table3_ephysSUA_task_stat_JCscript
% create and save a table ('Tephys') that contain new cells fields
% Fields:  Fr_BLE_Hz, SMA_MAX_Hz, SMA_MIN_Hz , Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,  ...
%      ttest_puf_H, ttest_del_H, ttest_res_H, ttest_puf_P, ttest_del_P,  ttest_res_P, ...
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Tephys')
% written by Julien Catanese 11/26/2018
% last updated JC 2/7/2019.
% last updated JC 3/22/2019: removed the Zscore part.
% last updated JC 4/12/2019: Changed name from Tephys to Tfig2_cor and Add Tcombo saving.

%% Get Firing Rate Values from FRepoch structure
Fr_BLE_Hz = FRepoch.BLE.mean;
Fr_puf_Hz = FRepoch.puf.mean;
Fr_del_Hz = FRepoch.del.mean;
Fr_res_Hz = FRepoch.res.mean;
SMA_MAX_Hz = max(SMA(:,200:end-200)')';
SMA_MIN_Hz = min(SMA(:,200:end-200)')';

%% Get ttest results from ttestEpoch structure
H_puf = ttestEpoch.H_puf_BLE;
H_del = ttestEpoch.H_del_BLE;
H_res = ttestEpoch.H_res_BLE;

P_puf = ttestEpoch.P_puf_BLE;
P_del = ttestEpoch.P_del_BLE;
P_res = ttestEpoch.P_res_BLE;

%% TTest ==> Define 7 Classes of cell response during task
tt_puf_Exc = (Tcombo.VMVL) &...
    (H_puf & (Fr_puf_Hz > Fr_BLE_Hz))   &...
    (~(H_del & (Fr_del_Hz < Fr_BLE_Hz))) &...
    (~(H_res & (Fr_res_Hz < Fr_BLE_Hz))) ;

tt_del_Exc = (Tcombo.VMVL) &...
    (H_del & (Fr_del_Hz > Fr_BLE_Hz))    &...
    (~(H_puf & (Fr_puf_Hz < Fr_BLE_Hz))) &...
    (~(H_res & (Fr_res_Hz < Fr_BLE_Hz))) ;

tt_res_Exc  = (Tcombo.VMVL) &...
    (H_res & (Fr_res_Hz > Fr_BLE_Hz))    &...
    (~(H_del & (Fr_del_Hz < Fr_BLE_Hz))) &...
    (~(H_puf & (Fr_puf_Hz < Fr_BLE_Hz))) ;

tt_type7_Exc = (Tcombo.VMVL) & tt_puf_Exc & tt_del_Exc & tt_res_Exc;
tt_type6_Exc = (Tcombo.VMVL) & tt_puf_Exc & tt_res_Exc & ~tt_type7_Exc;
tt_type5_Exc = (Tcombo.VMVL) & tt_del_Exc & tt_res_Exc & ~tt_type7_Exc;
tt_type4_Exc = (Tcombo.VMVL) & tt_puf_Exc & tt_del_Exc & ~tt_type7_Exc;
tt_type3_Exc = (Tcombo.VMVL) & tt_res_Exc & ~tt_del_Exc & ~tt_puf_Exc;
tt_type2_Exc = (Tcombo.VMVL) & tt_del_Exc & ~tt_res_Exc & ~tt_puf_Exc;
tt_type1_Exc = (Tcombo.VMVL) & tt_puf_Exc & ~tt_del_Exc & ~tt_res_Exc;
tt_type0     = (Tcombo.VMVL) & H_puf==0 &  H_del==0 & H_res==0;

%% Z-score ==> Define 7 Classes of cell response during task
Zthr = 3;
% z-score = (x-mean)/(std/sqrt(n))
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2); zSMA(find(zSMA==inf))=NaN;

[zbool zvalu] = pub_comp_Zthr_75msEpoch_JCfun(Zthr, zSMA, parfig)

zSMA_MAX_puf = zvalu.puf.max; zSMA_MAX_del = zvalu.del.max; zSMA_MAX_res = zvalu.res.max;
zSMA_MIN_puf = zvalu.puf.min; zSMA_MIN_del = zvalu.del.min; zSMA_MIN_res = zvalu.res.min;
zSMA_AbsPeak_all = max(abs(zSMA(:,200:end-200))')';
z_exct_pufAll = zbool.puf.exct; z_exct_delAll = zbool.del.exct; z_exct_resAll = zbool.res.exct;
z_inib_pufAll = zbool.puf.inib; z_inib_delAll = zbool.del.inib; z_inib_resAll = zbool.res.inib;

z_exctAll      =  z_exct_pufAll | z_exct_delAll | z_exct_resAll;   % All cell that are sig Excited in task
z_inibAll      =  z_inib_pufAll | z_inib_delAll | z_inib_resAll;
z_exct         =  z_exctAll & ~z_inibAll;  % Excite ONLY
z_inib         = ~z_exctAll &  z_inibAll; % Inib ONLY

z_complex      =  z_exctAll &  z_inibAll; % Complex ONLY
z_nosig        = ~z_exctAll & ~z_inibAll; % NoSig ONLY

z_exct_UniMod_1puf      =  z_exct &  z_exct_pufAll & ~z_exct_delAll  & ~z_exct_resAll;
z_exct_UniMod_2del      =  z_exct & ~z_exct_pufAll &  z_exct_delAll  & ~z_exct_resAll;
z_exct_UniMod_3res      =  z_exct & ~z_exct_pufAll & ~z_exct_delAll  &  z_exct_resAll;
z_inib_UniMod_1puf      =  z_inib &  z_inib_pufAll & ~z_inib_delAll  & ~z_inib_resAll;
z_inib_UniMod_2del      =  z_inib & ~z_inib_pufAll &  z_inib_delAll  & ~z_inib_resAll;
z_inib_UniMod_3res      =  z_inib & ~z_inib_pufAll & ~z_inib_delAll  &  z_inib_resAll;

z_exct_biMod_12      =  z_exct &  z_exct_pufAll &  z_exct_delAll  & ~z_exct_resAll;
z_exct_biMod_13      =  z_exct &  z_exct_pufAll & ~z_exct_delAll  &  z_exct_resAll;
z_exct_biMod_23      =  z_exct & ~z_exct_pufAll &  z_exct_delAll  &  z_exct_resAll;
z_inib_biMod_12      =  z_inib &  z_inib_pufAll &  z_inib_delAll  & ~z_inib_resAll;
z_inib_biMod_13      =  z_inib &  z_inib_pufAll & ~z_inib_delAll  &  z_inib_resAll;
z_inib_biMod_23      =  z_inib & ~z_inib_pufAll &  z_inib_delAll  &  z_inib_resAll;

z_exct_triMod_123      =  z_exct &  z_exct_pufAll &  z_exct_delAll  &  z_exct_resAll;
z_inib_triMod_123      =  z_inib &  z_inib_pufAll &  z_inib_delAll  &  z_inib_resAll;


%% SAVING TABLE
if parfig.saveTABLE ==1
    T1 = [];   T1=removevars(Tfig1_VMopto,[1,10,11,12,18,19,20,21,22 ])
    Tfig2_cor = addvars(T1, SMA_MAX_Hz, SMA_MIN_Hz,...
        Fr_BLE_Hz, Fr_puf_Hz, Fr_del_Hz, Fr_res_Hz,...
        zSMA_AbsPeak_all, zSMA_MAX_puf, zSMA_MAX_del,...
        zSMA_MAX_res, zSMA_MIN_puf, zSMA_MIN_del, zSMA_MIN_res,...
        z_exct_pufAll, z_exct_delAll, z_exct_resAll, ...
        z_inib_pufAll, z_inib_delAll, z_inib_resAll, ...
        z_exctAll, z_inibAll, z_exct, z_inib, z_complex, z_nosig, ...
        z_exct_UniMod_1puf, z_exct_UniMod_2del, z_exct_UniMod_3res, ...
        z_inib_UniMod_1puf, z_inib_UniMod_2del, z_inib_UniMod_3res, ...
        z_exct_biMod_12, z_exct_biMod_13, z_exct_biMod_23, ...
        z_inib_biMod_12, z_inib_biMod_13, z_inib_biMod_23, ...
        z_exct_triMod_123, z_inib_triMod_123,...
        H_puf, H_del, H_res, P_puf, P_del,  P_res,...
        tt_puf_Exc, tt_del_Exc, tt_res_Exc, tt_type0, tt_type1_Exc, tt_type2_Exc, ...
        tt_type3_Exc, tt_type4_Exc, tt_type5_Exc, tt_type6_Exc, tt_type7_Exc);
    
    save('listcell.mat','Tfig2_cor', '-append')
    disp('Tfig2_cor SAVED'); Tfig2_cor(1,:)
    
    % TCOMBO
%     Tcombo=[];
%     Tcombo=listcell(:,1:5); AreaID=Tcombo.AreaID; VMVL=Tcombo.VMVL;   
%     ncell =  Tcombo.ncell; nSess =  Tcombo.nSess; nMouse = Tcombo.nMouse;
%     Tcombo = addvars(Tcombo, nMouse, nSess ,VMVL, AreaID);    
%     Opto_inib = Topto.Opto_inib; Opto_exct = Topto.Opto_exct; Opto_post_sess= Topto.Opto_post_sess;
%     Tcombo= addvars(Tcombo, Opto_inib, Opto_exct, Opto_post_sess);   
%     z_exct=Tfig2_cor.z_exct; z_inib=Tfig2_cor.z_inib; z_complex=Tfig2_cor.z_complex;
%     z_nosig=Tfig2_cor.z_nosig; z_exct_UniMod_1puf=Tfig2_cor.z_exct_UniMod_1puf
%     z_exct_UniMod_2del=Tfig2_cor.z_exct_UniMod_2del; z_exct_UniMod_3res=Tfig2_cor.z_exct_UniMod_3res;
%     z_exct_biMod_12=Tfig2_cor.z_exct_biMod_12; z_exct_biMod_13=Tfig2_cor.z_exct_biMod_13;
%     z_exct_biMod_23=Tfig2_cor.z_exct_biMod_23; z_exct_triMod_123=Tfig2_cor.z_exct_triMod_123;    
%     Tcombo= addvars(Tcombo, z_exct, z_inib, z_complex, z_nosig, ...
%         z_exct_UniMod_1puf, z_exct_UniMod_2del, z_exct_UniMod_3res, ...
%         z_exct_biMod_12, z_exct_biMod_13, z_exct_biMod_23, ...
%         z_exct_triMod_123);
%     save([mypath '\Tcombo.mat'],'Tcombo'); Tcombo(1,:)
%     disp('Tcombo SAVED')
    
end



%% plot VENN DIAGRAM
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