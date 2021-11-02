function pub_fig_RasterSDF_zscoretypes_JCfun(Tcombo_z, SDF, listcell, parfig)
% function pub_fig_RasterSDF_zscoretypes_JCfun(Tcombo, SDF, parfig)
% plot Raster and Spike density function for single units (SUA). 
% Screen using Zscore for example of each type. 
% Written by Julien Catanese 3/15/2019
% last updated 4/1/19 JC
 
%% Id bentify each Type (using Ttest (FrBLE vs FrEpoch))  

% NEED to create TdiffLR using model from Tcombo_z

type0     = Tcombo_z.VMVL & TdiffLR_z.z_cont_nosig; 
type1_con = Tcombo_z.VMVL & TdiffLR_z.z_cont_UniMod_1puf;  
type2_con = Tcombo_z.VMVL & TdiffLR_z.z_cont_UniMod_2del;  
type3_con = Tcombo_z.VMVL & TdiffLR_z.z_cont_UniMod_3res;   
type4_con = Tcombo_z.VMVL & TdiffLR_z.z_cont_biMod_12pd;
type5_con = Tcombo_z.VMVL & TdiffLR_z.z_cont_biMod_23dr;
type6_con = Tcombo_z.VMVL & TdiffLR_z.z_cont_biMod_13pr;
type7_con = Tcombo_z.VMVL & TdiffLR_z.z_cont_triMod_123; 

Ntype0 = sum(type0)
Ntype1_con = sum(type1_con)
Ntype2_con = sum(type2_con)
Ntype3_con = sum(type3_con)
Ntype4_con = sum(type4_con)
Ntype5_con = sum(type5_con)
Ntype6_con = sum(type6_con)
Ntype7_con = sum(type7_con)

%% Type1 PUFF 
% Screening
parfig.title=['Type1 (#' num2str(Ntype1_con) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo_z(type1_con,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type1_con, SDF, parfig)
% Example  for final Figure 
listc1 = listcell('vgat15w10d3S1Ch4clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);


%% Type2 DELAY
% Screening
parfig.title=['Type2 (#' num2str(Ntype2_con) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo_z(type2_con,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type2_con, SDF, parfig)
% Example for final Figure 
listc1= listcell('vgat17w10d4S3Ch4clu#01', :);
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%% Type 3 RESP 
% Screening
parfig.title=['Type3 (#' num2str(Ntype3_con) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo_z(type3_con,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type3_con, SDF, parfig)
% Example for final Figure 
listc1= listcell('vgat11w10d4S4Ch6clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%% Type 4 (puff+delay)
% Screening
parfig.title=['Type4 (#' num2str(Ntype4_con) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo_z(type4_con,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type4_con, SDF, parfig)
% Example for final Figure  
listc1= listcell('vgat15w10d7S1Ch2clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%% Type 5  (RAMP)
% Screening
parfig.title=['Type5 (#' num2str(Ntype5_con) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo_z(type5_con,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type5_con, SDF, parfig)
% Example for final Figure  
listc1= listcell('vgat12w11d5S3Ch6clu#02', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

% listc1=[]; listc1= Tcombo_z('vgat12w11d5S3Ch2clu#02', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
%% Type 6  (BI-MOD)
% Screening
parfig.title=['Type6 (#' num2str(Ntype6_con) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo_z(type6_con,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type6_con, SDF, parfig)

listc1= listcell('vgat15w10d3S4Ch5clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%% Type 7 (TRI-MOD)
% Screening
parfig.title=['Type7 (#' num2str(Ntype7_con) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo_z(type7_con,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type7_con, SDF, parfig)

listc1= listcell('vgat15w10d7S1Ch1clu#02', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

listc1= listcell('vgat17w10d7S2Ch6clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%% Type 0  
% Screening
parfig.title=['Type0 (#' num2str(Ntype0) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo_z(type0,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type0, SDF, parfig)
%%
listc1= listcell('vgat14w14d5S2Ch6clu#01', :)
% listc1= listcell('vgat11w10d4S3Ch5clu#03',:)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

% vgat15w10d7S1Ch8clu#02
% vgat14w14d5S2Ch6clu#01

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Examples PuffCell with RASTER for final Figure
% listc1=[]; listc1= Tcombo_z('vgat15w10d7S1Ch2clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
% 
% listc1=[];listc1= Tcombo_z('vgat17w10d8S4Ch6clu#03', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
% 
% listc1=[];listc1= Tcombo_z('vgat15w10d3S1Ch4clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%%%%%%%%%%%%%%%%%%%%%%%%%
% % Examples DelayCell with RASTER for final Figure
% listc1= Tcombo_z('vgat17w10d4S3Ch4clu#01', :);
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Examples REspCell with RASTER for final Figure
% listc1= Tcombo_z('vgat11w10d4S4Ch6clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
% %
% listc1= Tcombo_z('vgat12w11d5S3Ch6clu#02', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
% %
% listc1= Tcombo_z('vgat17w10d8S3Ch1clu#03', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%%%%%%%%%%%%%%%%%%%%%%%%
% BOTH CELLS TYPE 6
% listc1= Tcombo_z('vgat17w10d7S3Ch3clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
% 


