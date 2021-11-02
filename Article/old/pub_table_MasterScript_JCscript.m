% pub_table_MasterScript_JCscript
% by JC 2/13/2019

%%
clear all, close all,
mypath = 'D:\JC_Analysis'
% mypath = 'C:\Users\Julien\Documents\WORK\JC_Analysis'
cd(mypath)

%% Get all events and task trials
% pub_Behavior_MasterScript_JCscript
% cd(mypath);

%% Get the Table with the list of recorded cells
FigSave_ON = 0; 
% pub_table1_listcell_JCscript

%% Make Tcoord : Coordinates and VM_VL Attribution (semi-Auto)
FigSave_ON=0;
% pub_table2_Tcoord_JCscript; close all;

%% Make Topto to identify cell that respond to OPTO-POST
% pub_table3_Topto_JCscript; close all;

%% Make Tephys (correct trials) to define fr resp types based on zscore(Tephys_z) or ttest(Tephys_tt) 
clearvars -except mypath parfig
load('listcell.mat')
parfig.plot=0;
parfig.saveSMA=1;
parfig.saveTABLE =1;

parfig.center_evt =  'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.BaselineEpoch= [150:2150];
parfig.pre = parfig.BaselineEpoch(end) + 1500 % define how much time before zero 
parfig.post = 1500 % define how much time after zero (zero will be defined by trigtimes.cor)

parfig.trial_type = {'cor'};
parfig.col = {'k'};
% pub_comp_ephySMA_task_ttest_JCfun(listcell, parfig)
%%
% pub_table4_Tfig2_cor_JCscript; close all;

%% Make TcLvR (ipsi-contra trials) to define diff zRight-zLeft (dzSMA TcLvR_z)  
% parfig.trial_type = {'cCL'};
% parfig.col = {'r'};
% [SMA_CL FRepoch_CL ttestEpoch_CL] = pub_comp_ephySMA_task_ttest_JCfun(listcell, parfig);
% parfig.trial_type = {'cCR'};
% parfig.col = {'b'};
% [SMA_CR FRepoch_CR ttestEpoch_CR] = pub_comp_ephySMA_task_ttest_JCfun(listcell, parfig);
%%
% pub_table5_Tfig3_ipsicontra_JCscript; close all;

%% Make Timpomi (impulse-omissions trials) to define diff zImp-zCor (dzSMA TcLvR_z)  
parfig.trial_type = {'imp'};
parfig.col = {'m'};
[SMA_imp FRepoch_imp ttestEpoch_imp] = pub_comp_ephySMA_task_ttest_JCfun(listcell, parfig);
parfig.trial_type =  {'omi'};
parfig.col = {'c'};
[SMA_omi FRepoch_omi ttestEpoch_omi] = pub_comp_ephySMA_task_ttest_JCfun(listcell, parfig);
%%
pub_table6_Timpomi_taskSUAimpomi_stat_JCscript; close all

% % Z-SCORE CLASSIFICATION : create Tephys_z and TcLvR_z (Z-Trh > 3std)
% pub_subT_TaskMod_Zscore_JCscript
% %
% pub_table6_Tcombo_JCscript; close all;

%% Figure1: Plot Behavior + Histo + Opto MUA
% pub_fig1_Behavior_JCscript;
% cd('D:\JC_Analysis');
% FigSave_ON=1
% pub_fig1_AllenAtlas_ShChan_3D_JCscript; close all;
% pub_disp_BrainAreaOpto_JCscript; close all;

%% Figure2: Ephys during corr trials (SUA and pop)
close all, clear all,
load('listcell.mat'); load('Tcombo'); 
load('SMA_cor_GoCue.mat');
parfig=parfigUsed;
parfig.center_evt =  'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.trial_type = {'cor'};
parfig.col = {'k'};
parfig.xlim = [-2750 1250];

zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2); zSMA(find(zSMA==inf))=NaN;

%% Fig 2A Venn Diagrams (2 types: +/-/complex/nosig and 7celltypes)
%colormap ( %'jet'   %  'flag' % 'copper'  % 'cool' % 'hot' % 'colorcube'   % 'bone'  % 'autumn' )
colmap = 'bone'
resolution = 0.05
pub_fig2_venndiagrams_zcor_JCfun(Tfig2_cor, colmap, resolution)


%% Fig 2B-I : Examples SUA => SDF-RASTER for each Type (from I to VII)
close all
parfig.saveSMA=0;
parfig.saveTABLE=0;
parfig.plot=1;
parfig.plotmerge=0;
parfig.ylabel = 'Fr (Hz)'
parfig.plotshaded= 'sem'
parfig.k=2;
pub_fig2_example_RasterSDF_JCfun(Tcombo, SMA, listcell, parfig)
% pub_fig_RasterSDF_ttesttypes_JCfun(SMA, parfig)

%% Fig 2J opto pie
pub_fig2_OptoPie_JCscript

%% Fig 2K-N : scatter Opto
pub_fig2_scatterOpto_JCscript

%% Fig 2O Distribution Boostrap Opto
parfig.Zthr = 3 ;
iter = 100;
parfig.title = 'Dist bootstrap'; parfig.ylim=[0 90];
pub_fig2_Distrib_TaskMod_Bootstrap_JCfun(zSMA, Tcombo, iter, parfig)

%% Fig 2P-Q  Matrice Populations
zSMA_VMexct= zSMA(logical(Tcombo.VMVL) & Tcombo.z_exct,:);
ZorTT= 'z_thr'
parfig.colormap = 'jet';
parfig.caxis = [-1 1];
parfig.sort_Xepoch = [parfig.pre-1500 : parfig.pre+1000];
parfig.trial_type = parfig.trial_type;
sort_variable = 'peak';
pub_fig2_MatricePop_JCfun(ZorTT, zSMA_VMexct, Tcombo, parfig)
% saveas()
parfig.sort_variable = 'Z-tr'
parfig.caxis = [0 9];
parfig.sort_Xepoch = [parfig.pre-750 : parfig.pre+1000];
pub_fig2_MatricePop_JCfun(ZorTT, zSMA_VMexct, Tcombo, parfig)
% saveas()

%% Fig 2R : Population firing rate
close all,
VMVL_all = Tcombo.VMVL;
Exct_only = Tcombo.z_exct & Tcombo.VMVL;
Inib_only = Tcombo.z_inib  & Tcombo.VMVL;
complex_only = Tcombo.z_complex  & Tcombo.VMVL;
nosig_only = Tcombo.z_nosig  & Tcombo.VMVL;
idx_list = logical([VMVL_all, Exct_only, Inib_only, complex_only, nosig_only ]);
parfig.ylabel = 'Zscore'
parfig.title = 'grand average'
parfig.col = {'k','r','g','c','y'}
parfig.plotshaded='sem'
parfig.k=2;
% parfig.plotshaded='std'
parfig.xlim = [0-2250 0+750];
parfig.ylim= [-5 7];
pub_fig_SMA_MUA_1by1_JCfun(idx_list, zSMA, parfig)
pub_fig_SMA_MUA_JCfun(logical(idx_list), zSMA, parfig)
legend ('VMVL all', ['(n=' num2str(sum(VMVL_all)) ')'],...
    'Exct only', ['(n=' num2str(sum(Exct_only)) ')'],...
    'inib only',  ['(n=' num2str(sum(Inib_only)) ')'], ...
    'complex', ['(n=' num2str(sum(complex_only)) ')'],...
    'nosig', ['(n=' num2str(sum(nosig_only)) ')'],...
    'Location', 'Best')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3 EPYS IPSI CONTRA (opto-Post)
clear all, close all,
load('listcell.mat');
load('SMA_cor_GoCue.mat');
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2); zSMA(find(zSMA==inf))=NaN;

parfig= parfigUsed;
parfig.saveSMA=0;
parfig.saveTABLE=0;
parfig.center_evt =  'GoCue';
parfig.xlim = [-2750 1250];
parfig.merge=1;
for iside = 1:2
    % IPSI
    if iside == 1
        load('SMA_cCL_GoCue.mat');
        parfigipsi.trial_type = {'cCL'};
        parfigipsi.col = {'r'};
        zSMAipsi = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
        SMAipsi= SMA; SSemAipsi=SSemA; clear SMA SSemA SSA;
        Tipsi = Tcomboipsi;
        % CONTRA
    elseif iside == 2
        load('SMA_cCR_GoCue.mat');
        parfigcont.trial_type = {'cCR'};
        parfigcont.col = {'b'};
        zSMAcont = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
        SMAcont= SMA; SSemAcont=SSemA; clear SMA SSemA SSA;
        Tcontra = Tcombocont;
    end
end

Zthr=3;
dzSMA = zSMAcont-zSMAipsi;
[zbool zvalu] = pub_subT_TaskMod_ZmaxThrEpoch_JCfun(Zthr, dzSMA, parfig)
% SAVE Sub_TABLE
T1 = removevars(Tcombo,'Var22');
Tdz_LR = addvars(T1(:,1:12), zbool.puf.exct, zbool.del.exct, zbool.res.exct,...
    zbool.puf.inib, zbool.del.inib, zbool.res.inib, zvalu.puf.max, ...
    zvalu.del.max, zvalu.res.max, zvalu.puf.min, zvalu.del.min, zvalu.res.min, ...
    'NewVariableNames', {'cont_puf' 'cont_del', 'cont_res','ipsi_puf', ...
    'ipsi_del', 'ipsi_res', 'diffZ_MAX_puf','diffZ_MAX_del','diffZ_MAX_res',...
    'diffZ_MIN_puf','diffZ_MIN_del','diffZ_MIN_res'}); Tdz_LR(1,:)
save('D:\JC_Analysis\listcell.mat','Tdz_LR', '-append'); disp('Tdz_LR SAVED');

contra_dz = Tcombo.z_exct & Tcombo.VMVL & (Tfig3_contipsi_dz.dzRL_exct_pufAll | Tdz_LR.cont_del | Tdz_LR.cont_res) ;
ipsi_dz = Tcombo.z_exct & Tcombo.VMVL & (Tdz_LR.ipsi_puf | Tdz_LR.ipsi_del | Tdz_LR.ipsi_res);

%% Fig 3A Venn Diagrams (2 types: +/-/complex/nosig and 7celltypes)
close all, colmap = 'jet';  resolution = 0.05;
% using zscore
Ncomplex_dz = sum(contra_dz & ipsi_dz)
data= [sum(contra_dz)-Ncomplex_dz , Ncomplex_dz, sum(ipsi_dz) - Ncomplex_dz]
vennX(data, resolution, colmap)
colormapeditor

% using ttest
% Ncontra_tt = sum(Tcombo.z_exct & Tcombo.VMVL & Tcombo.contra_cell & ~Tcombo.ipsi_cell);
% Nipsi_tt = sum(Tcombo.z_exct & Tcombo.VMVL & ~Tcombo.contra_cell & Tcombo.ipsi_cell);
% Ncomplex_tt = sum(Tcombo.z_exct & Tcombo.VMVL & Tcombo.contra_cell & Tcombo.ipsi_cell);
% data= [Ncontra_tt , Ncomplex_tt, Nipsi_tt]
% vennX(data, resolution, colmap)

%% Fig 3B Bar Graphs
% all cells
fig_bar_RvL_7types_JCscript % count ipsi only and contra only but not complex
fig_bar_RvL_3epochDZ_opto_JCscript
fig_bar_RvL_3epochDZ_Nosig_JCscript

close all, figure,
% contra_dz = Tcombo.z_exct & Tcombo.VMVL & (Tdz_LR.cont_puf | Tdz_LR.cont_del | Tdz_LR.cont_res) ;
data_mean =[abs(mean(zvalu.puf.max(contra_dz))), abs(mean(zvalu.puf.min(ipsi_dz))); ...
    abs(mean(zvalu.del.max(contra_dz))), abs(mean(zvalu.del.min(ipsi_dz)));...
    abs(mean(zvalu.res.max(contra_dz))), abs(mean(zvalu.res.min(ipsi_dz)))]
data_std = [abs(std(zvalu.puf.max(contra_dz))), abs(std(zvalu.puf.min(ipsi_dz))); ...
    abs(std(zvalu.del.max(contra_dz))), abs(std(zvalu.del.min(ipsi_dz)));...
    abs(std(zvalu.res.max(contra_dz))), abs(std(zvalu.res.min(ipsi_dz)))]
data_sem = [abs(std(zvalu.puf.max(contra_dz))/sqrt(size(contra_dz,1))), abs(std(zvalu.puf.min(ipsi_dz))/sqrt(size(ipsi_dz,1))); ...
    abs(std(zvalu.del.max(contra_dz))/sqrt(size(contra_dz,1))), abs(std(zvalu.del.min(ipsi_dz))/sqrt(size(ipsi_dz,1)));...
    abs(std(zvalu.res.max(contra_dz))/sqrt(size(contra_dz,1))), abs(std(zvalu.res.min(ipsi_dz))/sqrt(size(ipsi_dz,1)))]

ranksum(zvalu.puf.max(contra_dz), zvalu.puf.min(ipsi_dz))
bar(data_mean)
hold on,
errorbar([1-0.15 1+0.15; 2-0.15 2+0.15; 3-0.15 3+0.15] ,data_mean, data_sem*3, '.')

%% Fig 3C-H : Examples SUA => SDF-RASTER for each Type (from I to VII)
close all
parfig.saveSMA=0;
parfig.saveTABLE=0;
parfig.plotmerge=1;
parfig.nRaster = 2;
parfig.xlim= [-1750 1000];
parfig.plot=1;
parfig.k=2;
parfig.ylabel = 'z';
parfig.plotshaded= 'sem';
parfig.parfigcont=parfigcont;
parfig.parfigipsi = parfigipsi;
dzSMA=  zSMAcont-zSMAipsi;
pub_fig3_example_RasterSDF_ipsicontra_JCfun(Tcombo, dzSMA, listcell, parfig)
% saveas(gcf, ['D:\JC_Figures\SDF_RASTER_cor\SDF_RASTER_cor_listcell_ncell#' num2str(ncell) ],'jpg')

%% Fig 3I-J Mean DZ
Xt=[1:1:size(dzSMA,2)]-parfig.pre;
parfig.xlim=[-2000 750]
K=2;

figure,
hold on, plot(Xt, mean(dzSMA(contra_dz,:)))
hold on, plot(Xt, mean(dzSMA(ipsi_dz,:)))

pub_fig3_LvR_shaded_dzSMA_JCscript

%% Fig 3K-L Distribution Bootstrap Opto
% contra (Fig 3K)
parfig.Zthr = 2;
parfig.bootstrap.iter = 500
parfig.bootstrap.groupsize=40
parfig.title = 'CONTRA Dist bootstrap'; parfig.ylim=[0 70];
group_bool= [contra_dz & Tcombo.Opto_post_sess &  ~(Tcombo.Opto_inib | Tcombo.Opto_exct), ...
    contra_dz & Tcombo.Opto_post_sess & Tcombo.Opto_inib, ...
    contra_dz & Tcombo.Opto_post_sess & Tcombo.Opto_exct ];
group_ID = {'conf int limit (3std)';...
    ['Thal NOL cells (' num2str(sum(Tcombo.VMVL & Tcombo.Opto_post_sess & Tcombo.z_exct & ~(Tcombo.Opto_inib | Tcombo.Opto_exct))) ')'];...
    ['SvTh NEG cells (0' num2str(num2str(sum(Tcombo.VMVL & Tcombo.Opto_post_sess & Tcombo.z_exct & Tcombo.Opto_inib))) ')' ];...
    ['SvTh POS cells (0' num2str(num2str(sum(Tcombo.VMVL & Tcombo.Opto_post_sess & Tcombo.z_exct & Tcombo.Opto_exct))) ')' ]};
pub_fig3_Distrib_Bootstrap_JCfun(dzSMA, group_ID, group_bool, parfig)
xlim(parfig.xlim);

% ipsi (Fig 3L)
parfig.Zthr = 2;
parfig.bootstrap.iter = 500
parfig.bootstrap.groupsize=40
parfig.title = ' IPSI Dist bootstrap'; parfig.ylim=[0 70];
group_bool= [ipsi_dz & Tcombo.Opto_post_sess &  ~(Tcombo.Opto_inib | Tcombo.Opto_exct), ...
    ipsi_dz & Tcombo.Opto_post_sess & Tcombo.Opto_inib, ...
    ipsi_dz & Tcombo.Opto_post_sess & Tcombo.Opto_exct ];
GroupID = {'conf int limit (3std)';...
    ['Thal NOL cells (' num2str(sum(Tcombo.VMVL & Tcombo.Opto_post_sess & Tcombo.z_exct & ~(Tcombo.Opto_inib | Tcombo.Opto_exct))) ')'];...
    ['SvTh NEG cells (0' num2str(num2str(sum(Tcombo.VMVL & Tcombo.Opto_post_sess & Tcombo.z_exct & Tcombo.Opto_inib))) ')' ];...
    ['SvTh POS cells (0' num2str(num2str(sum(Tcombo.VMVL & Tcombo.Opto_post_sess & Tcombo.z_exct & Tcombo.Opto_exct))) ')' ]};
pub_fig3_Distrib_Bootstrap_JCfun(-dzSMA, group_ID, group_bool, parfig)
xlim(parfig.xlim);

%% Fig 3M-N Matrices ipsi vs Contra
close all,
dzSMA = zSMAcont-zSMAipsi;
dzSMA_contra = abs(dzSMA(contra_dz,:));
dzSMA_ipsi = abs(dzSMA(ipsi_dz,:));

ZorTT= 'z_thr'
parfig.sort_variable = 'peak';
parfig.colormap = 'jet';
parfig.caxis = [-1 1];
parfig.sort_Xepoch = [parfig.pre-1500 : parfig.pre+750];
parfig.trial_type = parfigcont.trial_type;
pub_fig2_MatricePop_JCfun(ZorTT, -dzSMA_contra, Tcombo, parfig)
parfig.trial_type = parfigipsi.trial_type;
pub_fig2_MatricePop_JCfun(ZorTT, dzSMA_ipsi, Tcombo, parfig)

parfig.sort_variable = 'Z-tr'
parfig.colormap = 'jet';
parfig.caxis = [-4.5 4.5];
parfig.sort_Xepoch = [parfig.pre-1500 : parfig.pre+1000];
parfig.trial_type = parfigcont.trial_type;
pub_fig2_MatricePop_JCfun(ZorTT, -dzSMA_contra, Tcombo, parfig)
parfig.trial_type = parfigipsi.trial_type;
pub_fig2_MatricePop_JCfun(ZorTT, dzSMA_ipsi, Tcombo, parfig)
pub_fig2_MatricePop_JCfun(ZorTT, -dzSMA(Tcombo.z_exct & Tcombo.VMVL,:), Tcombo, parfig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fig4: Plot OPTO Behavior
ddisp('Fig4 to do')
% pub_fig_BehavOpto_BarMeanStat_JCscript

%% Fig5: Plot Ephys Impusle Omission Opto Example

clear all, close all,
load('listcell.mat');
load('SMA_cor_GoCue.mat');
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2); zSMA(find(zSMA==inf))=NaN;

parfig= parfigUsed;
parfig.saveSMA=0;
parfig.saveTABLE=0;
parfig.center_evt =  'GoCue';
parfig.xlim = [-2750 1250];
parfig.merge=1;
for iside = 1:2
    % IPSI
    if iside == 1
        load('SMA_cCL_GoCue.mat');
        parfigipsi.trial_type = {'cCL'};
        parfigipsi.col = {'r'};
        zSMAipsi = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
        SMAipsi= SMA; SSemAipsi=SSemA; clear SMA SSemA SSA;
        Tipsi = Tcomboipsi;
        % CONTRA
    elseif iside == 2
        load('SMA_cCR_GoCue.mat');
        parfigcont.trial_type = {'cCR'};
        parfigcont.col = {'b'};
        zSMAcont = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
        SMAcont= SMA; SSemAcont=SSemA; clear SMA SSemA SSA;
        Tcontra = Tcombocont;
    end
end

Zthr=3;
dzSMA = zSMAcont-zSMAipsi;
[zbool zvalu] = pub_subT_TaskMod_ZmaxThrEpoch_JCfun(Zthr, dzSMA, parfig)
% SAVE Sub_TABLE
T1 = removevars(Tcombo,'Var22');
Tdz_LR = addvars(T1(:,1:12), zbool.puf.exct, zbool.del.exct, zbool.res.exct,...
    zbool.puf.inib, zbool.del.inib, zbool.res.inib, zvalu.puf.max, ...
    zvalu.del.max, zvalu.res.max, zvalu.puf.min, zvalu.del.min, zvalu.res.min, ...
    'NewVariableNames', {'cont_puf' 'cont_del', 'cont_res','ipsi_puf', ...
    'ipsi_del', 'ipsi_res', 'diffZ_MAX_puf','diffZ_MAX_del','diffZ_MAX_res',...
    'diffZ_MIN_puf','diffZ_MIN_del','diffZ_MIN_res'}); Tdz_LR(1,:)
save('D:\JC_Analysis\listcell.mat','Tdz_LR', '-append'); disp('Tdz_LR SAVED');

contra_dz = Tcombo.z_exct & Tcombo.VMVL & (Tdz_LR.cont_puf | Tdz_LR.cont_del | Tdz_LR.cont_res) ;
ipsi_dz = Tcombo.z_exct & Tcombo.VMVL & (Tdz_LR.ipsi_puf | Tdz_LR.ipsi_del | Tdz_LR.ipsi_res);
















%% Fig6: Plot Ephys opto task
% Get SMA for opto trial
parfig.trial_type = {'opt'};
parfig.saveTABLE = 0;
parfig.saveSMA=1;
parfig.plot=0;

pub_fig6_EphysOptoTask_JCscript

%% plot Classifier
