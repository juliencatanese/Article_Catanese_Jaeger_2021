% Classifier_Average_Left_Right_JCscript
%  written by Julien Catanese 10/31/2018
% last updated: 10/31/2018

%%
clear all, close all,

MouseID = 'JCVGAT1';
FolderID = dir(['D:\JC_Analysis\' MouseID '*\*taskopto*']);

ALLM1=[]; ALLM2=[]; ALLM3=[];

for nf=2:max(size(FolderID));
    disp([num2str(nf) '/' num2str(max(size(FolderID))) '   MouseID =' MouseID  '   FolderID = ' FolderID(nf).name]);
    cd([FolderID(nf).folder '\' FolderID(nf).name]);

    
load('ClassM1_APuff_cCL_cCR.mat'); AP1=M1; M1=[];  
load('ClassM1_Delay_cCL_cCR.mat'); DL1=M1; M1=[];  
load('ClassM1_GoCue_cCL_cCR.mat'); GC1=M1; M1=[];  

if max(size(AP1))>40
    AP1=AP1(1:40)
    DL1=DL1(1:40)
    GC1=GC1(1:40)
end
V=ones(1,40-max(size(AP1)))*nan ;

VAP1=[AP1 V];
VDL1=[DL1 V];
VGC1=[GC1 V];

ALLM1 = [ALLM1; VAP1];
ALLM2 = [ALLM2; VDL1];
ALLM3 = [ALLM3; VGC1];

end

ALLM1_mean = nanmean(ALLM1)
ALLM2_mean = nanmean(ALLM2)
ALLM3_mean = nanmean(ALLM3)

ALLM1_std = nanstd(ALLM1)
ALLM2_std = nanstd(ALLM2)
ALLM3_std = nanstd(ALLM3)

figure, hold on, 
plot(ALLM1_mean, 'r')
plot(ALLM2_mean, 'b')
plot(ALLM3_mean, 'g')
% plot(ALLM1_mean+ALLM1_std, '--m')
% plot(ALLM1_mean-ALLM1_std, '--m')
% plot(ALLM2_mean+ALLM2_std, '--c')
% plot(ALLM2_mean-ALLM2_std, '--c')
% plot(ALLM3_mean+ALLM3_std, '--g')
% plot(ALLM3_mean-ALLM3_std, '--g')

% plot(ALLM3_mean, 'g')