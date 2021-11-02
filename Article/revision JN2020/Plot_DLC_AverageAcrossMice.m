% Plot_DLC_AverageAcrossMice  
% Julien Catanese 10/19/2020

%% Single trial Analysis 
clc, close all, clear all, 
% Mouse#3 Sess#1
FileLocation{1}= 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d8_z4900_VM_taskopto_optopost_G912_180715_vidY_150tr_20cel_13mW_tooDeep';
MouseID =  'vgat14w14d8';
cd(FileLocation{1}); 
Get_Abs_Average_3DLCvar_JCScript
plot_Abs_AverageLR_3DLCvar_JCScript
plot_Abs_Average_3DLCvar_JCScript2

% Mouse#4 Sess#2
MouseID =  'vgat15w10d7';
FileLocation{2}= 'D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d7_z4400_VM_taskopto_optopost_G912_180721_vidY_250tr_24cel_11mW';
cd(FileLocation{2});
Get_Abs_Average_3DLCvar_JCScript
plot_Abs_AverageLR_3DLCvar_JCScript
plot_Abs_Average_3DLCvar_JCScript2

% Mouse#5 Sess#3
FileLocation{3} = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW';
MouseID =  'vgat17w10d7';
cd(FileLocation{3}); 
Get_Abs_Average_3DLCvar_JCScript
plot_Abs_AverageLR_3DLCvar_JCScript
plot_Abs_Average_3DLCvar_JCScript2

% Mouse#5 Sess#4
MouseID =  'vgat17w10d5';
FileLocation{4}= 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d5_z4300_VM_taskopto_optopost_EOBD_180728_vidY_200tr_20cel_12mW_Q4';
cd(FileLocation{4}); 
Get_Abs_Average_3DLCvar_JCScript
plot_Abs_AverageLR_3DLCvar_JCScript
plot_Abs_Average_3DLCvar_JCScript2


%% Average across animal 
clc
for iSess = 1:4
   cd(FileLocation{iSess}); 
   load('Average_ABS_DLC.mat')
   MV1(iSess,:) = MeanABS1; MV2(iSess,:) = MeanABS2; MV3(iSess,:) = MeanABS3;
   SV1(iSess,:) = SemABS1 ; SV2(iSess,:) = SemABS2 ; SV3(iSess,:) = SemABS3 ;
   
end
figure, hold on, 
plot(mean(MV1),'g','LineWidth',2)
plot(mean(MV2),'c','LineWidth',2)
plot(mean(MV3),'y','LineWidth',2)
xline(1.5*25,'k--','LineWidth',2);
K=1
X=[1:1:100];
plotshaded(X, [mean(MV1)-mean(SV1,1); mean(MV1)+mean(SV1,1)],'g');
plotshaded(X, [mean(MV2)-mean(SV2,1); mean(MV2)+mean(SV2,1)],'c');
plotshaded(X, [mean(MV3)-mean(SV3,1); mean(MV3)+mean(SV3,1)],'y');