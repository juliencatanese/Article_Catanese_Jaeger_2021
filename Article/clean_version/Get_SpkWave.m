% Get_SpkWaves
% Julien Catanese 11/5/2020

close all 


%% Fig2 cell#156 RAMP VM
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d7_z4350_VM_taskonly_nonepost_CC4F_180714_vidN_075tr_40cel_00mW_pre_drift';
cd(FileLocation)
CellNum=156;
load('times_S1Ch4_sub.mat') 
clu=2 
C1= find(cluster_class(:,1)==clu);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')




%% Fig2 cell#109     % vgat14w14d2S2Ch1clu#01
CellNum = 109;
cd('D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d2_z4670_VM_taskopto_optopost_G912_180709_vidY_100tr_29cel_10mW');
load('times_S2Ch1_sub.mat') 
C1= find(cluster_class(:,1)==1);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')


%% Fig2 cell#284
CellNum = 284;
cd('D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d7_z4400_VM_taskopto_optopost_G912_180721_vidY_250tr_24cel_11mW');
load('times_S1Ch2_sub.mat') 
C1= find(cluster_class(:,1)==1);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')

%% Fig2 cell#39 
cd('D:\DATA EMORY\JC_Analysis\JCVGAT11\vgat11_w10d4_z4300_VM_taskonly_nonepost_CC4F_180505_vidY_150tr_40cel_00mW')
load('times_S4Ch6_sub.mat') 
CellNum= 39
C1= find(cluster_class(:,1)==1);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')

%% Fig2 cell#76 
cd('D:\DATA EMORY\JC_Analysis\JCVGAT12\vgat12_w11d5_z4300_VM_taskopto_optopost_CCAE_180609_vidY_100tr_42cel_10mW_3otr')
load('times_S3Ch6_sub.mat') 
CellNum= 76
C1= find(cluster_class(:,1)==2);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')

%% Fig2 cell#243
cd('D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d3_z4250_VM_taskonly_nonepost_CAED_180717_vidY_150tr_37cel_00mW')
load('times_S1Ch4_sub.mat') 
CellNum= 243
C1= find(cluster_class(:,1)==1);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')


%% Fig2 cell# 278
cd('D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d3_z4250_VM_taskonly_nonepost_CAED_180717_vidY_150tr_37cel_00mW')
load('times_S4Ch5_sub.mat') 
CellNum= 278
C1= find(cluster_class(:,1)==1);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')


%% Fig2 cell# 281
cd('D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d7_z4400_VM_taskopto_optopost_G912_180721_vidY_250tr_24cel_11mW')
load('times_S1Ch1_sub.mat') 
CellNum= 281
C1= find(cluster_class(:,1)==2);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')

%% Fig2 cell#398
cd('D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d4_z4200_VM_taskopto_optopost_CCAE_180726_vidY_075tr_41cel_10mW_bl_4otr')
load('times_S3Ch4_sub.mat') 
CellNum= 398
C1= find(cluster_class(:,1)==1);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')

%% Fig4 cell#74 
cd('D:\DATA EMORY\JC_Analysis\JCVGAT12\vgat12_w11d5_z4300_VM_taskopto_optopost_CCAE_180609_vidY_100tr_42cel_10mW_3otr')
load('times_S3Ch5_sub.mat') 
CellNum= 74
C1= find(cluster_class(:,1)==2);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')

%% Fig4 cell#198
cd('D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d8_z4900_VM_taskopto_optopost_G912_180715_vidY_150tr_20cel_13mW_tooDeep')
load('times_S3Ch2_sub.mat') 
CellNum= 198
C1= find(cluster_class(:,1)==1);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')

%% Fig4 cell#377
cd('D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d4_z4200_VM_taskopto_optopost_CCAE_180726_vidY_075tr_41cel_10mW_bl_4otr')
load('times_S2Ch3_sub.mat') 
CellNum=377 
clu=1 
C1= find(cluster_class(:,1)==clu);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')

%% Fig2 cell#15
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT11\vgat11_w10d4_z4300_VM_taskonly_nonepost_CC4F_180505_vidY_150tr_40cel_00mW';
cd(FileLocation)
CellNum=15;
load('times_S2Ch7_sub.mat') 
clu=1 
C1= find(cluster_class(:,1)==clu);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')

%% Fig2 cell#273
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d3_z4250_VM_taskonly_nonepost_CAED_180717_vidY_150tr_37cel_00mW';
cd(FileLocation)
CellNum=273;
load('times_S4Ch4_sub.mat') 
clu=1 
C1= find(cluster_class(:,1)==clu);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')



%% Fig2 cell#156 RAMP
FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d7_z4350_VM_taskonly_nonepost_CC4F_180714_vidN_075tr_40cel_00mW_pre_drift';
cd(FileLocation)
CellNum=156;
load('times_S1Ch4_sub.mat') 
clu=2 
C1= find(cluster_class(:,1)==clu);
K=1;
AV = mean(spikes(C1,:))
SD = std(spikes(C1,:))
figure,  plot(mean(spikes(C1,:)),'k')
hold on,  plotshaded([1:1:40], [AV-(SD); AV+(SD)] , 'k')
ylim([-150 50])
title(['cell#' num2str(CellNum) ' (shaded=1STD)'])
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'emf')
saveas(gcf, ['C:\Users\catan\Dropbox\Article_JNeuro_Julien_Dieter\REVISION\ISI and rpv\SpkWaves_Cell#' num2str(CellNum)],'png')


