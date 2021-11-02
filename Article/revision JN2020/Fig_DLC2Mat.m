close all, clear all 
cd ('D:\DATA EMORY\JC_Analysis')
load Table_DLC_vgat14w12d2.mat

%%
close all
figure, hold on, 
subplot(3,1,1), plot(allTab.NoseX); hold on; plot(allTab.NoseY); title('NoseX')
subplot(3,1,2), plot(allTab.MouthX); hold on; plot(allTab.MouthY); title('MouthX')
subplot(3,1,3),  plot(allTab.TongueProba,'r'); title('TongueProba')

figure, hold on, 
plot(allTab.NoseY, 'b'); 
plot(allTab.MouthY, 'g')
plot(allTab.TongueProba*100,'r')
legend('NoseY','MouthY','TongueProba')

figure, hold on, 
plot(allTab.TongueProba*100,'r')
plot(allTab.TongueY,'c')
legend('Tongue proba', 'TongueY' )

%% average over trials (100 frames)
NoseX_tr = allTab.NoseX(1:100)
for itr = 1:(size(allTab,1)/100)-1
    NoseX_tr = [NoseX_tr allTab.NoseX(itr*100:((itr+1)*100)-1)];
end

NoseY_tr = allTab.NoseY(1:100)
for itr = 1:(size(allTab,1)/100)-1
    NoseY_tr = [NoseY_tr allTab.NoseY(itr*100:((itr+1)*100)-1)];
end

MouthX_tr = allTab.MouthX(1:100)
for itr = 1:(size(allTab,1)/100)-1
    MouthX_tr = [MouthX_tr allTab.MouthX(itr*100:((itr+1)*100)-1)];
end

MouthY_tr = allTab.MouthY(1:100)
for itr = 1:(size(allTab,1)/100)-1
    MouthY_tr = [MouthY_tr allTab.MouthY(itr*100:((itr+1)*100)-1)];
end

TongueX_tr = allTab.TongueX(1:100)
for itr = 1:(size(allTab,1)/100)-1
    TongueX_tr = [TongueX_tr allTab.TongueX(itr*100:((itr+1)*100)-1)];
end

TongueY_tr = allTab.TongueY(1:100)
for itr = 1:(size(allTab,1)/100)-1
    TongueY_tr = [TongueY_tr allTab.TongueY(itr*100:((itr+1)*100)-1)];
end

TongueProba_tr = allTab.TongueProba(1:100)
for itr = 1:(size(allTab,1)/100)-1
    TongueProba_tr = [TongueProba_tr allTab.TongueProba(itr*100:((itr+1)*100)-1)];
end


%% plot averaged trials and save  
close all 
VAR_DLC = TongueY_tr
VARname_str = 'TongueY'
Fig_Average_LvR_VAR_DLC(MouthY_tr, TongueProba_tr, VAR_DLC, VARname_str)

VAR_DLC = NoseY_tr
VARname_str = 'NoseY'
Fig_Average_LvR_VAR_DLC(MouthY_tr, TongueProba_tr, VAR_DLC, VARname_str)


