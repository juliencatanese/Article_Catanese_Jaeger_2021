% pub_fig2_scatterOpto_JCscript
% scatter plot firing rate for each type of cells (in the Excited pool) 
% compare each task epochs to baseline  
% color indicate the cell opto modulated (blue=inhibited, magenta=excited) 
% Written by Julien Catanese 3/21/2019

close all 
puf_exct = Tcombo.VMVL & (Tcombo_z.z_exct_UniMod_1puf | Tcombo_z.z_exct_biMod_12pd);
del_exct = Tcombo.VMVL & (Tcombo_z.z_exct_UniMod_2del | Tcombo_z.z_exct_biMod_23dr) ;
res_exct = Tcombo.VMVL & (Tcombo_z.z_exct_UniMod_3res | Tcombo_z.z_exct_biMod_13pr) ;
tri_exct = Tcombo.VMVL & Tcombo_z.z_exct_triMod_123; 
all_exct = Tcombo.VMVL & Tcombo_z.z_exct & Tcombo.Opto_post_sess

% PIE OPTO (for figure inset)
figure, explode=[0 0 1]
pie([sum(puf_exct & Tcombo.Opto_inib),sum(puf_exct & Tcombo.Opto_exct), sum(puf_exct & Tcombo.Opto_post_sess & ~(Tcombo.Opto_exct | Tcombo.Opto_inib))], explode)
title('PUFF CELL')

figure, explode=[0 0 1]
pie([sum(del_exct & Tcombo.Opto_inib),sum(del_exct & Tcombo.Opto_exct), sum(del_exct & Tcombo.Opto_post_sess & ~(Tcombo.Opto_exct | Tcombo.Opto_inib))], explode)
title('DELAY CELL')

figure, explode=[0 0 1]
pie([sum(res_exct & Tcombo.Opto_inib),sum(res_exct & Tcombo.Opto_exct), sum(res_exct & Tcombo.Opto_post_sess & ~(Tcombo.Opto_exct | Tcombo.Opto_inib))], explode)
title('RESPONSE CELL')

figure, explode=[0 0 1]
pie([sum(tri_exct & Tcombo.Opto_inib), sum(tri_exct & Tcombo.Opto_exct), sum(tri_exct & Tcombo.Opto_post_sess & ~(Tcombo.Opto_exct | Tcombo.Opto_inib))], explode)
title('TRI CELL')

%% SCATTER OPTO 
close all, 
parfig.title = 'SAMPLE CELL'
parfig.XlimMAX = 20;
figure , 
scatter(Tephys.Fr_BLE_Hz(all_exct), Tephys.Fr_puf_Hz(all_exct), 60, [0.8 0.8 0.8], 'diamond', 'filled') 
hold on, 
scatter(Tephys.Fr_BLE_Hz(puf_exct), Tephys.Fr_puf_Hz(puf_exct), 60, 'k', 'diamond', 'filled'), 
hold on , 
scatter(Tephys.Fr_BLE_Hz(puf_exct & Tcombo.Opto_exct), Tephys.Fr_puf_Hz(puf_exct & Tcombo.Opto_exct), 60, 'm', 'diamond', 'filled'), 
hold on, 
scatter(Tephys.Fr_BLE_Hz(puf_exct & Tcombo.Opto_inib), Tephys.Fr_puf_Hz(puf_exct & Tcombo.Opto_inib), 60, 'c',  'diamond', 'filled'), 
hold on, 
plot([0:1:100],[0:1:100],'k' ,'LineWidth', 1),
xlabel('baseline (Hz)'), 
ylabel('sample (Hz)') 
xlim([0 parfig.XlimMAX]), 
ylim([0 parfig.XlimMAX]); 
title(parfig.title); 
pause(0.1)
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatterOpto_PUFCELL_pufpEpoch', 'png')
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatterOpto_PUFCELL_pufpEpoch', 'emf')
pause(0.1)

%

parfig.title = 'DELAY CELL'
parfig.XlimMAX = 20;
figure , 
% scatter(Tephys.Fr_BLE_Hz(logical(Tcombo.VMVL)), Tephys.Fr_del_Hz(logical(Tcombo.VMVL)), 60, [0.8 0.8 0.8], 'filled') 
% hold on, 
scatter(Tephys.Fr_BLE_Hz(del_exct), Tephys.Fr_del_Hz(del_exct), 60, 'k','filled'), 
hold on , 
scatter(Tephys.Fr_BLE_Hz(del_exct & Tcombo.Opto_exct), Tephys.Fr_del_Hz(del_exct & Tcombo.Opto_exct), 60, 'm', 'filled'), 
hold on, 
scatter(Tephys.Fr_BLE_Hz(del_exct & Tcombo.Opto_inib), Tephys.Fr_del_Hz(del_exct & Tcombo.Opto_inib), 60, 'c', 'filled'), 
hold on, 
plot([0:1:100],[0:1:100],'k' ,'LineWidth', 1),
xlabel('baseline (Hz)'), 
ylabel('delay (Hz)') 
xlim([0 parfig.XlimMAX]), 
ylim([0 parfig.XlimMAX]); 
title(parfig.title); 
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatterOpto_DELCELL_delpEpoch', 'png')
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatterOpto_DELCELL_delpEpoch', 'emf')
pause(0.1)


%
parfig.title = 'RESPONSE CELL'
parfig.XlimMAX = 20;
figure , 
% scatter(Tephys.Fr_BLE_Hz(logical(Tcombo.VMVL)), Tephys.Fr_res_Hz(logical(Tcombo.VMVL)), 60, [0.8 0.8 0.8], 'filled') 
% hold on, 
scatter(Tephys.Fr_BLE_Hz(res_exct), Tephys.Fr_res_Hz(res_exct), 60, 'k','filled'), 
hold on , 
scatter(Tephys.Fr_BLE_Hz(res_exct & Tcombo.Opto_exct), Tephys.Fr_res_Hz(res_exct & Tcombo.Opto_exct), 60, 'm', 'filled'), 
hold on, 
scatter(Tephys.Fr_BLE_Hz(res_exct & Tcombo.Opto_inib), Tephys.Fr_res_Hz(res_exct & Tcombo.Opto_inib), 60, 'c', 'filled'), 
hold on, 
plot([0:1:100],[0:1:100],'k' ,'LineWidth', 1),
xlabel('baseline (Hz)'), 
ylabel('response (Hz)') 
xlim([0 parfig.XlimMAX]), 
ylim([0 parfig.XlimMAX]); 
title(parfig.title); 
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatterOpto_RESCELL_respEpoch', 'png')
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatterOpto_RESCELL_respEpoch', 'emf')
pause(0.1)

parfig.title = 'TRI CELL'
parfig.XlimMAX = 15;
figure , 
FR_TASK = mean([Tephys.Fr_del_Hz Tephys.Fr_puf_Hz Tephys.Fr_res_Hz],2); 
% scatter(Tephys.Fr_BLE_Hz(logical(Tcombo.VMVL)), FR_TASK(logical(Tcombo.VMVL)), 60, [0.8 0.8 0.8], 'filled') 
% hold on, 
scatter(Tephys.Fr_BLE_Hz(all_exct), FR_TASK(all_exct), 60, 'k','filled'),
hold on , 
scatter(Tephys.Fr_BLE_Hz(all_exct & Tcombo.Opto_exct), FR_TASK(all_exct & Tcombo.Opto_exct), 60, 'm', 'filled'), 
hold on, 
scatter(Tephys.Fr_BLE_Hz(all_exct & Tcombo.Opto_inib), FR_TASK(all_exct & Tcombo.Opto_inib), 60, 'c', 'filled'), 
hold on, 
plot([0:1:100],[0:1:100],'k' ,'LineWidth', 1),
xlabel(' Baseline (Hz)'), 
ylabel('Task Average (Hz)') 
xlim([0 parfig.XlimMAX]), 
ylim([0 parfig.XlimMAX]); 
title(parfig.title); 
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatterOpto_ALLCELL_TaskEpoch', 'png')
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatterOpto_ALLCELL_Task  Epoch', 'emf')
pause(0.1)

