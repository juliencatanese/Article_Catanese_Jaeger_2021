% pub_fig5_example_impomi_Raster1by1_JCscript 
% Plot rasters for 3 example cells with IMP and OMI 
% listc1= listcell('vgat12w11d5S3Ch5clu#02', :)  
% listc1= listcell('vgat14w14d8S3Ch2clu#01', :)  
% listc1= listcell('vgat17w10d4S2Ch3clu#01', :) 
% Each raster ploted separately will have to be recombined in Illustrator
% Written by Julien Catanese 4/16/2019 

close all; 

Folder2Save = 'D:\JC_Figures\ARTICLE_JC_DJ\Fig5\'

% listc1= listcell('vgat12w11d5S3Ch5clu#02', :)  
% listc1= listcell('vgat14w14d8S3Ch2clu#01', :)  
% listc1= listcell('vgat17w10d4S2Ch3clu#01', :)  

% CELL#1 from VGAT12
listc1= listcell('vgat12w11d5S3Ch5clu#02', :)  % cell 74 
figure, 
parfig.trial_type = {'imp'}; 
parfig.col = {'m'}
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
ylim([0 20]); 
saveas(gcf, [Folder2Save 'Fig5C_EXimp_vgat12w11d5S3Ch5clu#02'], 'emf')
saveas(gcf, [Folder2Save 'Fig5C_EXimp_vgat12w11d5S3Ch5clu#02'], 'png')
figure, 
parfig.trial_type = {'omi'}; 
parfig.col = {'c'}
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
ylim([0 20]); 
saveas(gcf, [Folder2Save 'Fig5C_EXomi_vgat12w11d5S3Ch5clu#02'], 'emf')
saveas(gcf, [Folder2Save 'Fig5C_EXomi_vgat12w11d5S3Ch5clu#02'], 'png')
figure, 
parfig.trial_type = {'cor'}; 
parfig.col = {'k'}
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
ylim([0 20]); 
saveas(gcf, [Folder2Save 'Fig5C_EXcor_vgat12w11d5S3Ch5clu#02'], 'emf')
saveas(gcf, [Folder2Save 'Fig5C_EXcor_vgat12w11d5S3Ch5clu#02'], 'png')


% CELL#2 from VGAT14 
listc1= listcell('vgat14w14d8S3Ch2clu#01', :)     % cell 198
figure, 
parfig.trial_type = {'imp'}; 
parfig.col = {'m'}
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
ylim([0 15]); 
saveas(gcf, [Folder2Save 'Fig5D_EXimp_vgat14w14d8S3Ch2clu#01'], 'emf')
saveas(gcf, [Folder2Save 'Fig5D_EXimp_vgat14w14d8S3Ch2clu#01'], 'png')
figure, 
parfig.trial_type = {'omi'}; 
parfig.col = {'c'}
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
ylim([0 15]); 
saveas(gcf, [Folder2Save 'Fig5D_EXomi_vgat14w14d8S3Ch2clu#01'], 'emf')
saveas(gcf, [Folder2Save 'Fig5D_EXomi_vgat14w14d8S3Ch2clu#01'], 'png')
figure, 
parfig.trial_type = {'cor'}; 
parfig.col = {'k'}
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
ylim([0 15]); 
saveas(gcf, [Folder2Save 'Fig5D_EXcor_vgat14w14d8S3Ch2clu#01'], 'emf')
saveas(gcf, [Folder2Save 'Fig5D_EXcor_vgat14w14d8S3Ch2clu#01'], 'png')


% CELL#3 from VGAT17 
listc1= listcell('vgat17w10d4S2Ch3clu#01', :)  % cell 377
figure, 
parfig.trial_type = {'imp'}; 
parfig.col = {'m'}
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
ylim([0 25]); 
saveas(gcf, [Folder2Save 'Fig5E_EXimp_vgat17w10d4S2Ch3clu#01'], 'emf')
saveas(gcf, [Folder2Save 'Fig5E_EXimp_vgat17w10d4S2Ch3clu#01'], 'png')
figure, 
parfig.trial_type = {'omi'}; 
parfig.col = {'c'}
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
ylim([0 25]); 
saveas(gcf, [Folder2Save 'Fig5E_EXomi_vgat17w10d4S2Ch3clu#01'], 'emf')
saveas(gcf, [Folder2Save 'Fig5E_EXomi_vgat17w10d4S2Ch3clu#01'], 'png')
figure, 
parfig.trial_type = {'cor'}; 
parfig.col = {'k'}
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
ylim([0 25]); 
saveas(gcf, [Folder2Save 'Fig5E_EXcor_vgat17w10d4S2Ch3clu#01'], 'emf')
saveas(gcf, [Folder2Save 'Fig5E_EXcor_vgat17w10d4S2Ch3clu#01'], 'png')





