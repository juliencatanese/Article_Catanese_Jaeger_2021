% PLOT_RASTER_PSTH_JCscript
%% GoCue: Raster SDF Shaded
Sh=str2num(CellID(2));Ch=str2num(CellID(5)); CLUST=str2num(CellID(9));
psth_trial_type={'cor'}; col={'k'};pre=1000; post=1000; K=1;  
%GaussSmooth=40; 
psth_center_evt='GoCue'; 
[sdf_mean, sdf_sem, nbtrial] = pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'png')
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'emf')




%% Licks: Raster SDF Shaded
% psth_center_evt='Licks'; 
% pub_SDF_Raster_1cell_ShadedSEM_JCfun(Sh, Ch, CLUST, psth_center_evt , psth_trial_type, K, col, pre, post, GaussSmooth)
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'png')
% saveas(gcf, ['Average_SDF_Example_CellSh' num2str(Sh) 'Ch' num2str(Ch) 'clu'  num2str(CLUST) '_'  psth_center_evt],'emf')
