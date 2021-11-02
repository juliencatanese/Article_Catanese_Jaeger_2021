% Get_DLC_VAR_Trial_Mat_script
%% average over trials (100 frames)
NoseX_tr = allTab.NoseX(1:100);
for itr = 1:(size(allTab,1)/100)-1;
    NoseX_tr = [NoseX_tr allTab.NoseX(itr*100:((itr+1)*100)-1)];
end

NoseY_tr = allTab.NoseY(1:100);
for itr = 1:(size(allTab,1)/100)-1;
    NoseY_tr = [NoseY_tr allTab.NoseY(itr*100:((itr+1)*100)-1)];
end

MouthX_tr = allTab.MouthX(1:100);
for itr = 1:(size(allTab,1)/100)-1;
    MouthX_tr = [MouthX_tr allTab.MouthX(itr*100:((itr+1)*100)-1)];
end

MouthY_tr = allTab.MouthY(1:100);
for itr = 1:(size(allTab,1)/100)-1;
    MouthY_tr = [MouthY_tr allTab.MouthY(itr*100:((itr+1)*100)-1)];
end

TongueX_tr = allTab.TongueX(1:100);
for itr = 1:(size(allTab,1)/100)-1;
    TongueX_tr = [TongueX_tr allTab.TongueX(itr*100:((itr+1)*100)-1)];
end

TongueY_tr = allTab.TongueY(1:100);
for itr = 1:(size(allTab,1)/100)-1;
    TongueY_tr = [TongueY_tr allTab.TongueY(itr*100:((itr+1)*100)-1)];
end

TongueProba_tr = allTab.TongueProba(1:100);
for itr = 1:(size(allTab,1)/100)-1;
    TongueProba_tr = [TongueProba_tr allTab.TongueProba(itr*100:((itr+1)*100)-1)];
end

WhiskerX_tr = allTab.WhiskerX(1:100);
for itr = 1:(size(allTab,1)/100)-1;
    WhiskerX_tr = [WhiskerX_tr allTab.WhiskerX(itr*100:((itr+1)*100)-1)];
end

WhiskerY_tr = allTab.WhiskerY(1:100);
for itr = 1:(size(allTab,1)/100)-1;
    WhiskerY_tr = [WhiskerY_tr allTab.WhiskerY(itr*100:((itr+1)*100)-1)];
end

WhiskerProba_tr = allTab.WhiskerProba(1:100);
for itr = 1:(size(allTab,1)/100)-1;
    WhiskerProba_tr = [WhiskerProba_tr allTab.WhiskerProba(itr*100:((itr+1)*100)-1)];
end

display('DLC var sorted by trials: completed')
