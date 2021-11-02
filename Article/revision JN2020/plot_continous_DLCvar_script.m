% plot_continous_DLCvar_script  

figure, hold on, title('X vs Y')
subplot(3,1,1), plot(allTab.NoseX); hold on; plot(allTab.NoseY); title('Nose X vs Y')
subplot(3,1,2), plot(allTab.MouthX); hold on; plot(allTab.MouthY); title('Mouth X vs Y')
subplot(3,1,3), plot(allTab.TongueX); hold on; plot(allTab.TongueY); title('Tongue X vs Y')


figure, hold on, 
plot(allTab.NoseY, 'b'); 
plot(allTab.MouthY, 'g')
plot(allTab.TongueProba*100,'r')
legend('NoseY','MouthY','TongueProba')
title('Timing all Y vs Tongue proba')

figure, hold on, 
plot(allTab.TongueY,'c')
plot(allTab.TongueProba*100,'r')
legend('TongueY','Tongue proba'  )
title('TONGUE Y vs Proba')

figure, hold on, 
plot(allTab.TongueY, 'r')
plot(allTab.NoseY, 'b'); 
plot(allTab.MouthY, 'g')
legend('TongueY','NoseY','MouthY')
title('all Y VAR')

figure, hold on, 
plot(allTab.NoseY, 'b'); 
plot(allTab.NoseX*3,'c')
plot(allTab.MouthY, 'g')
plot(allTab.MouthX, 'y')
legend('NoseY','NoseX','MouthY','MouthX')
title('mouth nose XY coVAR') 

% plot all trial Average
figure, hold on,
plot(mean(TongueProba_tr,2)*100), plot(mean(TongueY_tr,2)), plot(mean(TongueX_tr,2))
legend('Tongue proba', 'TongueY', 'TongueX')
title('Tongue all trial average DLC')

figure, hold on,
plot(mean(MouthY_tr,2)), plot(mean(MouthX_tr,2)), plot(mean(NoseY_tr,2)), plot(mean(NoseX_tr,2))
legend('MouthY', 'MouthX', 'NoseY','NoseX')
title('Nose and Mouth all trial average DLC')

figure, hold on,
plot(mean(WhiskerY_tr,2)), plot(mean(WhiskerX_tr,2)), plot(mean(WhiskerProba_tr,2)*100)
legend('WhiskerY', 'WhiskerX', 'Whisker Proba')
title('Whisker all trial average DLC')

%% OPTIONAL:  plot Left vs Rigth trials average and save
% idxL_dlc: Select L vs R trials based on MouthY
% VAR_DLC=MouthY_tr;
% VAR_str = 'MouthY_tr'
% Thr = median(mean(VAR_DLC));
% AV=[];
% for itr=1:125;
%     idxT = find(TongueProba_tr(:,itr)>0.95);
%     AV = [AV mean(VAR_DLC(idxT,itr))];
% end
% idxL_dlc = find(AV>Thr);
% idxR_dlc = find(AV<Thr);
% 
% idxL_dlc = find(AV>Thr);
% idxR_dlc = find(AV<Thr);
% idx_str1 = 'Ldlc'
% idx_str2 = 'Rdlc'
% 
% plot_Average_trialtype_DLCvar(VAR_DLC, idxL_dlc, idxR_dlc, VAR_str, idx_str1, idx_str2)
