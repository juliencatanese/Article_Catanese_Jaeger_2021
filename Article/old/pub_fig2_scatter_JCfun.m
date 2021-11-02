function pub_fig2_scatter_JCfun(Tephys_type, parfig)
% function pub_fig2_scatter_JCfun(Tephys_type)
% INPUT= Table Tephys, possibility to select for a given group of cells.
% OUTPUT= 3 scatters plot for Figure2 E,F,G,H pub_fig2_scatter_JCfun  
% Written by Julien Catanese 03/13/2019


XLIMmax = parfig.XlimMAX; 

figure , scatter(Tephys_type.Fr_BLE_Hz, Tephys_type.Fr_puf_Hz,  10, 'k'), 
xlabel('baseline (Hz)'), ylabel('sample (Hz)') 
hold on, plot([0:1:100],[0:1:100],'r' ,'LineWidth', 2),
xlim([0 XLIMmax]), ylim([0 XLIMmax]); title(parfig.title); 

figure , scatter(Tephys_type.Fr_BLE_Hz, Tephys_type.Fr_del_Hz, 10, 'k'), 
xlabel('baseline (Hz)'), ylabel('delay (Hz)') 
hold on, plot([0:1:100],[0:1:100],'r' ,'LineWidth', 2),
xlim([0 XLIMmax]), ylim([0 XLIMmax]);title(parfig.title); 

figure , scatter(Tephys_type.Fr_BLE_Hz, Tephys_type.Fr_res_Hz, 10, 'k'), 
xlabel('baseline (Hz)'), ylabel('response (Hz)') 
hold on, plot([0:1:100],[0:1:100],'r' ,'LineWidth', 2),
xlim([0 XLIMmax]), ylim([0 XLIMmax]);title(parfig.title); 





