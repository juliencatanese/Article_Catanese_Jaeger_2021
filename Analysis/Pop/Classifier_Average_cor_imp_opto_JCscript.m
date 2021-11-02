% Classifier_Average_cor_imp_opto_JCscript
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
    
load('ClassM1_cor_imp.mat');
load('ClassM2_cor_oNO.mat');
load('ClassM3_oNO_imp.mat');
V=ones(1,100-max(size(M1)))*nan 
V1=[M1 V]
V2=[M2 V]
V3=[M3 V]

ALLM1 = [ALLM1; V1];
ALLM2 = [ALLM2; V2];
ALLM3 = [ALLM3; V3];

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
plot(ALLM1_mean+ALLM1_std, '--r')
plot(ALLM1_mean-ALLM1_std, '--r')

plot(ALLM3_mean, 'g')