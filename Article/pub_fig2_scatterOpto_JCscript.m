% pub_fig2_scatterOpto_JCscript
% scatter plot firing rate for each type of cells (in the Excited pool) 
% compare each task epochs to baseline  
% color indicate the cell opto modulated (blue=inhibited, magenta=excited) 
% Written by Julien Catanese 3/21/2019

close all 
puf_exct = Tcombo.VMVL & (Tcombo.z_exct_UniMod_1puf | Tcombo.z_exct_biMod_12);
del_exct = Tcombo.VMVL & (Tcombo.z_exct_UniMod_2del | Tcombo.z_exct_biMod_23) ;
res_exct = Tcombo.VMVL & (Tcombo.z_exct_UniMod_3res | Tcombo.z_exct_biMod_13) ;
tri_exct = Tcombo.VMVL & Tcombo.z_exct_triMod_123; 
all_exct = Tcombo.VMVL & Tcombo.z_exct & Tcombo.Opto_post_sess
all_inib = Tcombo.VMVL & Tcombo.z_inib & Tcombo.Opto_post_sess

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

% SCATTER OPTO 
close all, 
parfig.title = 'SAMPLE CELL'
parfig.YlimMAX = 12;
parfig.XlimMAX = 12;
figure , 
% scatter(Tfig2_cor.Fr_BLE_Hz(all_exct), Tfig2_cor.Fr_puf_Hz(all_exct), 60, [0.8 0.8 0.8], 'diamond', 'filled') 
% hold on, 
scatter(Tfig2_cor.Fr_BLE_Hz(puf_exct), Tfig2_cor.Fr_puf_Hz(puf_exct), 60, [0.45 0.45 0.45], 'diamond', 'filled'), 
hold on , 
scatter(Tfig2_cor.Fr_BLE_Hz(puf_exct & Tcombo.Opto_exct), Tfig2_cor.Fr_puf_Hz(puf_exct & Tcombo.Opto_exct), 60, 'm', 'diamond', 'filled'), 
hold on, 
scatter(Tfig2_cor.Fr_BLE_Hz(puf_exct & Tcombo.Opto_inib), Tfig2_cor.Fr_puf_Hz(puf_exct & Tcombo.Opto_inib), 60, 'c',  'diamond', 'filled'), 
hold on, 
plot([0:1:100],[0:1:100],'k' ,'LineWidth', 1),
xlabel('baseline (Hz)'), 
ylabel('sample (Hz)') 
xlim([0 parfig.XlimMAX]), 
ylim([0 parfig.YlimMAX]); 
title(parfig.title); pause(0.1)
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatOptDiam_type1cell_SAMPLEepoch', 'png');
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatOptDiam_type1cell_SAMPLEepoch', 'emf');pause(0.1)
%
%
parfig.title = 'DELAY CELL'
figure , 
% scatter(Tfig2_cor.Fr_BLE_Hz(all_exct), Tfig2_cor.Fr_puf_Hz(all_exct), 60, [0.8 0.8 0.8], 'diamond', 'filled') 
% hold on, 
scatter(Tfig2_cor.Fr_BLE_Hz(del_exct), Tfig2_cor.Fr_del_Hz(del_exct), 60, [0.45 0.45 0.45], 'diamond','filled'), 
hold on , 
scatter(Tfig2_cor.Fr_BLE_Hz(del_exct & Tcombo.Opto_exct), Tfig2_cor.Fr_del_Hz(del_exct & Tcombo.Opto_exct), 60, 'm', 'diamond', 'filled'), 
hold on, 
scatter(Tfig2_cor.Fr_BLE_Hz(del_exct & Tcombo.Opto_inib), Tfig2_cor.Fr_del_Hz(del_exct & Tcombo.Opto_inib), 60, 'c', 'diamond', 'filled'), 
hold on, 
plot([0:1:100],[0:1:100],'k' ,'LineWidth', 1),
xlabel('baseline (Hz)'), 
ylabel('delay (Hz)') 
xlim([0 parfig.XlimMAX]), 
ylim([0 parfig.YlimMAX]); 
title(parfig.title); pause(0.1)
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatOptDiam_type2cell_DELAYepoch', 'png')
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatOptDiam_type2cell_DELAYepoch', 'emf')
pause(0.1)


%
parfig.title = 'RESPONSE CELL'
figure , 
% scatter(Tfig2_cor.Fr_BLE_Hz(all_exct), Tfig2_cor.Fr_puf_Hz(all_exct), 60, [0.8 0.8 0.8], 'diamond', 'filled') 
% hold on, 
scatter(Tfig2_cor.Fr_BLE_Hz(res_exct), Tfig2_cor.Fr_res_Hz(res_exct), 60, [0.45 0.45 0.45], 'diamond','filled'), 
hold on , 
scatter(Tfig2_cor.Fr_BLE_Hz(res_exct & Tcombo.Opto_exct), Tfig2_cor.Fr_res_Hz(res_exct & Tcombo.Opto_exct), 60, 'm', 'diamond', 'filled'), 
hold on, 
scatter(Tfig2_cor.Fr_BLE_Hz(res_exct & Tcombo.Opto_inib), Tfig2_cor.Fr_res_Hz(res_exct & Tcombo.Opto_inib), 60, 'c', 'diamond', 'filled'), 
hold on, 
plot([0:1:100],[0:1:100],'k' ,'LineWidth', 1),
xlabel('baseline (Hz)'), 
ylabel('response (Hz)') 
xlim([0 parfig.XlimMAX]), 
ylim([0 parfig.YlimMAX]); 
title(parfig.title); pause(0.1);
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatOptDiam_type3cell_RESPepoch', 'png')
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatOptDiam_type3cell_RESPepoch', 'emf')
pause(0.1)

parfig.title = 'TRI CELL'
figure , 
FR_TASK = mean([Tfig2_cor.Fr_del_Hz Tfig2_cor.Fr_puf_Hz Tfig2_cor.Fr_res_Hz],2); 
% scatter(Tfig2_cor.Fr_BLE_Hz(all_exct), Tfig2_cor.Fr_puf_Hz(all_exct), 60, [0.8 0.8 0.8], 'diamond', 'filled') 
% hold on, 
scatter(Tfig2_cor.Fr_BLE_Hz(all_exct), FR_TASK(all_exct), 60, [0.45 0.45 0.45], 'diamond','filled'),
hold on , 
scatter(Tfig2_cor.Fr_BLE_Hz(all_exct & Tcombo.Opto_exct), FR_TASK(all_exct & Tcombo.Opto_exct), 60, 'm', 'diamond', 'filled'), 
hold on, 
scatter(Tfig2_cor.Fr_BLE_Hz(all_exct & Tcombo.Opto_inib), FR_TASK(all_exct & Tcombo.Opto_inib), 60, 'c', 'diamond', 'filled'), 
hold on, 
plot([0:1:100],[0:1:100],'k' ,'LineWidth', 1),
xlabel(' Baseline (Hz)'), 
ylabel('Task Average (Hz)') 
xlim([0 parfig.XlimMAX]), 
ylim([0 parfig.YlimMAX]); 
title(parfig.title); pause(0.1); 
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatOptDiam_type4Tricell_TASKepoch', 'png')
saveas(gcf, 'D:\JC_Figures\ARTICLE_JC_DJ\Fig2\Fig2_ScatOptDiam_type4Tricell_TASKepoch', 'emf'); pause(0.1)

