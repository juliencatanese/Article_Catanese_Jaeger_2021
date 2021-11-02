function pub_fig2_example_RasterSDF_JCfun(Tcombo, SDF, listcell, parfig)
% function pub_fig_RasterSDF_zscoretypes_JCfun(Tcombo, SDF, parfig)
% plot Raster and Spike density function for single units (SUA). 
% Screen using Zscore for example of each type. 
% Written by Julien Catanese 3/15/2019
% last updated 4/1/19 JC
 
%% Identify each Type (using zscore) 
    
type0     = Tcombo.VMVL & Tcombo.z_nosig; 
type1_Exc = Tcombo.VMVL & Tcombo.z_exct_UniMod_1puf;  
type2_Exc = Tcombo.VMVL & Tcombo.z_exct_UniMod_2del;  
type3_Exc = Tcombo.VMVL & Tcombo.z_exct_UniMod_3res;   
type4_Exc = Tcombo.VMVL & Tcombo.z_exct_biMod_12;
type5_Exc = Tcombo.VMVL & Tcombo.z_exct_biMod_23;
type6_Exc = Tcombo.VMVL & Tcombo.z_exct_biMod_13;
type7_Exc = Tcombo.VMVL & Tcombo.z_exct_triMod_123; 


Ntype0 = sum(type0)
Ntype1_Exc = sum(type1_Exc)
Ntype2_Exc = sum(type2_Exc)
Ntype3_Exc = sum(type3_Exc)
Ntype4_Exc = sum(type4_Exc)
Ntype5_Exc = sum(type5_Exc)
Ntype6_Exc = sum(type6_Exc)
Ntype7_Exc = sum(type7_Exc)

%% test 
% listc1=[]; listc1= listcell('vgat17w10d3S2Ch5clu#01', :)
% if parfig.plotmerge==1, 
% figure, 
% parfig.col = parfig.parfigcont.col 
% parfig.trial_type = parfig.parfigcont.trial_type
% parfig.currentpanelNb = 1; 
% end
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
% 
% parfig.col = parfig.parfigipsi.col 
% parfig.trial_type = parfig.parfigipsi.trial_type
% parfig.currentpanelNb = 2; 
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%% Type1 PUFF 
% Screening
parfig.title=['Type1 (#' num2str(Ntype1_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type1_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type1_Exc, SDF, parfig)
% Example  for final Figure 
listc1 = listcell('vgat15w10d3S1Ch4clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);


%% Type2 DELAY
% Screening
parfig.title=['Type2 (#' num2str(Ntype2_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type2_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type2_Exc, SDF, parfig)
% Example for final Figure 
listc1= listcell('vgat17w10d4S3Ch4clu#01', :);
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%% Type 3 RESP 
% Screening
parfig.title=['Type3 (#' num2str(Ntype3_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type3_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type3_Exc, SDF, parfig)
% Example for final Figure 
listc1= listcell('vgat11w10d4S4Ch6clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%% Type 4 (puff+delay)
% Screening
parfig.title=['Type4 (#' num2str(Ntype4_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type4_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type4_Exc, SDF, parfig)
% Example for final Figure  
listc1= listcell('vgat15w10d7S1Ch2clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%% Type 5  (RAMP)
% Screening
parfig.title=['Type5 (#' num2str(Ntype5_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type5_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type5_Exc, SDF, parfig)
% Example for final Figure  
listc1= listcell('vgat12w11d5S3Ch6clu#02', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

% listc1=[]; listc1= Tcombo('vgat12w11d5S3Ch2clu#02', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
%% Type 6  (BI-MOD)
% Screening
parfig.title=['Type6 (#' num2str(Ntype6_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type6_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type6_Exc, SDF, parfig)

listc1= listcell('vgat15w10d3S4Ch5clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%% Type 7 (TRI-MOD)
% Screening
parfig.title=['Type7 (#' num2str(Ntype7_Exc) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type7_Exc,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type7_Exc, SDF, parfig)

listc1= listcell('vgat15w10d7S1Ch1clu#02', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

listc1= listcell('vgat17w10d7S2Ch6clu#01', :)
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%% Type 0  
% Screening
parfig.title=['Type0 (#' num2str(Ntype0) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type0,:), SDF, parfig)
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
% listc1=[]; listc1= Tcombo('vgat15w10d7S1Ch2clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
% 
% listc1=[];listc1= Tcombo('vgat17w10d8S4Ch6clu#03', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
% 
% listc1=[];listc1= Tcombo('vgat15w10d3S1Ch4clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%%%%%%%%%%%%%%%%%%%%%%%%%
% % Examples DelayCell with RASTER for final Figure
% listc1= Tcombo('vgat17w10d4S3Ch4clu#01', :);
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Examples REspCell with RASTER for final Figure
% listc1= Tcombo('vgat11w10d4S4Ch6clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
% %
% listc1= Tcombo('vgat12w11d5S3Ch6clu#02', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
% %
% listc1= Tcombo('vgat17w10d8S3Ch1clu#03', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%%%%%%%%%%%%%%%%%%%%%%%%
% BOTH CELLS TYPE 6
% listc1= Tcombo('vgat17w10d7S3Ch3clu#01', :)
% pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
% 


