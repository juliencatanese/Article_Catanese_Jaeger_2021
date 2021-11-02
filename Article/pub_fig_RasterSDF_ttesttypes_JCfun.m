function pub_fig_RasterSDF_types_JCfun(SDF, parfig)
% plot Raster and Spike density function for single units (SUA). 
% Screen using Zscore for example of each type. 
% Written by Julien Catanese 3/15/2019
 
load('listcell.mat')

%% Identify each Type (using Ttest (FrBLE vs FrEpoch)) 
puf_Exc = Tcombo.VMVL &...
    (Tephys.H_puf & (Tephys.Fr_puf_Hz > Tephys.Fr_BLE_Hz))    &...
    (~(Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz))) ; 

del_Exc = Tcombo.VMVL &...
    (Tephys.H_del & (Tephys.Fr_del_Hz > Tephys.Fr_BLE_Hz))    &...
    (~(Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz))) ; 

res_Exc  = Tcombo.VMVL &...
    (Tephys.H_res & (Tephys.Fr_res_Hz > Tephys.Fr_BLE_Hz))    &...
    (~(Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz))) ; 


type7_Exc = Tcombo.VMVL & puf_Exc & del_Exc & res_Exc;
type6_Exc = Tcombo.VMVL & puf_Exc & res_Exc & ~type7_Exc;
type5_Exc = Tcombo.VMVL & del_Exc & res_Exc & ~type7_Exc;
type4_Exc = Tcombo.VMVL & puf_Exc & del_Exc & ~type7_Exc;
type3_Exc = Tcombo.VMVL & res_Exc & ~del_Exc & ~puf_Exc; 
type2_Exc = Tcombo.VMVL & del_Exc & ~res_Exc & ~puf_Exc; 
type1_Exc = Tcombo.VMVL & puf_Exc & ~del_Exc & ~res_Exc; 
type0     = Tcombo.VMVL & Tephys.H_puf==0 &  Tephys.H_del==0 & Tephys.H_res==0; 

Ntype0 = sum(type0)
Ntype1_Exc = sum(type1_Exc)
Ntype2_Exc = sum(type2_Exc)
Ntype3_Exc = sum(type3_Exc)
Ntype4_Exc = sum(type4_Exc)
Ntype5_Exc = sum(type5_Exc)
Ntype6_Exc = sum(type6_Exc)
Ntype7_Exc = sum(type7_Exc)

%% test 
% listcell4=[]; listcell4= listcell('vgat15w10d3S1Ch6clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);

%% Type1 PUFF 
% Screening
parfig.title=['Type1 (#' num2str(Ntype1_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(listcell(type1_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type1_Exc, SDF, parfig)
% Example  for final Figure 
listcell4=[];listcell4= listcell('vgat15w10d3S1Ch4clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);


%% Type2 DELAY
% Screening
parfig.title=['Type2 (#' num2str(Ntype2_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(listcell(type2_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type2_Exc, SDF, parfig)
% Example for final Figure 
listcell4= listcell('vgat17w10d4S3Ch4clu#01', :);
pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);

%% Type 3 RESP 
% Screening
parfig.title=['Type3 (#' num2str(Ntype3_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(listcell(type3_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type3_Exc, SDF, parfig)
% Example for final Figure 
listcell4= listcell('vgat11w10d4S4Ch6clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);

%% Type 4 (puff+delay)
% Screening
parfig.title=['Type4 (#' num2str(Ntype4_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(listcell(type4_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type4_Exc, SDF, parfig)
% Example for final Figure  
listcell4=[]; listcell4= listcell('vgat15w10d7S1Ch2clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);

%% Type 5  (RAMP)
% Screening
parfig.title=['Type5 (#' num2str(Ntype5_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(listcell(type5_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type5_Exc, SDF, parfig)
% Example for final Figure  
listcell4=[]; listcell4= listcell('vgat12w11d5S3Ch6clu#02', :)
pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);

% listcell4=[]; listcell4= listcell('vgat12w11d5S3Ch2clu#02', :)
% pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);
%% Type 6  (BI-MOD)
% Screening
parfig.title=['Type6 (#' num2str(Ntype6_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(listcell(type6_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type6_Exc, SDF, parfig)

listcell4=[]; listcell4= listcell('vgat15w10d3S4Ch5clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);

%% Type 7 (TRI-MOD)
% Screening
parfig.title=['Type7 (#' num2str(Ntype7_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(listcell(type7_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type7_Exc, SDF, parfig)

listcell4= listcell('vgat15w10d7S1Ch1clu#02', :)
pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);

listcell4= listcell('vgat17w10d7S2Ch6clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);

%% Type 0  
% Screening
parfig.title=['Type0 (#' num2str(Ntype0) ')' ]
pub_fig_SMA_SUA_JCfun2(listcell(type0,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type0, SDF, parfig)
%%
listcell4= listcell('vgat14w14d5S2Ch6clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);

% vgat15w10d7S1Ch8clu#02
% vgat14w14d5S2Ch6clu#01

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Examples PuffCell with RASTER for final Figure
% listcell4=[]; listcell4= listcell('vgat15w10d7S1Ch2clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);
% 
% listcell4=[];listcell4= listcell('vgat17w10d8S4Ch6clu#03', :)
% pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);
% 
% listcell4=[];listcell4= listcell('vgat15w10d3S1Ch4clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);

%%%%%%%%%%%%%%%%%%%%%%%%%
% % Examples DelayCell with RASTER for final Figure
% listcell4= listcell('vgat17w10d4S3Ch4clu#01', :);
% pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Examples REspCell with RASTER for final Figure
% listcell4= listcell('vgat11w10d4S4Ch6clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);
% %
% listcell4= listcell('vgat12w11d5S3Ch6clu#02', :)
% pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);
% %
% listcell4= listcell('vgat17w10d8S3Ch1clu#03', :)
% pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);

%%%%%%%%%%%%%%%%%%%%%%%%
% BOTH CELLS TYPE 6
% listcell4= listcell('vgat17w10d7S3Ch3clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listcell4, parfig);
% 


