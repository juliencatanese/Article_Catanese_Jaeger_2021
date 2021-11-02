% plot_Abs_Average_3DLCvar_JCScript
% function plot_Abs_Average_DLCvar_JCfun(VAR1, VAR2, VARname_str1, VARname_str2)
% Julien Catanese 9/23/2020

load([ MouseID '_DLCresults.mat'])
load('Ntrial_type.mat')
load ('evt.mat','evt_trial', 'evt_lick_L', 'evt_lick_R')
load ('time.mat');
Get_DLC_VAR_Trial_Mat_script

%% Plot VAR using idxL and idxR
VAR1 = TongueY_tr; VAR_str1 = 'TongueY'; 
VAR2 = NoseY_tr; VAR_str2 = 'NoseY'; 
VAR3 = WhiskerX_tr; VAR_str3 = 'WhiskerX'; 
Thr1 = mean(mean(VAR1(10:30,:)));
Thr2 = mean(mean(VAR2(10:30,:)));
Thr3 = mean(mean(VAR3(10:30,:)));

idx_str = 'cor';
idx = sort([trial.idx_correct_L trial.idx_correct_R])
ntr=size(idx,2);

figure, hold on, 
K=1

MeanABS1 = mean(abs(VAR1(:,idx)-Thr1),2); 
MeanABS2 = mean(abs(VAR2(:,idx)-Thr2),2); 
MeanABS3 = mean(abs(VAR3(:,idx)-Thr3),2); 

SemABS1 = (std(abs(VAR1(:,idx)-Thr1)',1)/sqrt(ntr))'; 
SemABS2 = (std(abs(VAR2(:,idx)-Thr2)',1)/sqrt(ntr))'; 
SemABS3 = (std(abs(VAR3(:,idx)-Thr3)',1)/sqrt(ntr))'; 


X=[1:1:100]';
plot(MeanABS1,'g','LineWidth',2);
plot(MeanABS2,'c','LineWidth',2);
plot(MeanABS3,'y','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
plotshaded(X', [MeanABS1-(SemABS1*K) MeanABS1+(SemABS1*K)],'g');
plotshaded(X', [MeanABS2-(SemABS2*K) MeanABS2+(SemABS2*K)],'c');
plotshaded(X', [MeanABS3-(SemABS3*K) MeanABS3+(SemABS3*K)],'y');

legend(VAR_str1 ,VAR_str2, VAR_str3  ,'GoCue');
title( ['DLC events' MouseID]);

saveas(gcf, ['D:\DATA EMORY\JC_Analysis\Average_3DLCresults_LR' MouseID '_' VAR_str1 '_' VAR_str2 '_' VAR_str3 ],'png')
saveas(gcf, ['D:\DATA EMORY\JC_Analysis\Average_3DLCresults_LR' MouseID '_' VAR_str1 '_' VAR_str2 '_' VAR_str3 ],'emf')

save('Average_ABS_DLC.mat','VAR_str1','VAR_str2','VAR_str3',...
    'MeanABS1','MeanABS2','MeanABS3',...
    'SemABS1','SemABS2','SemABS3',...
    'Thr1','Thr2','Thr3');

