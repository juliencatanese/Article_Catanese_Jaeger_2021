% plot_Abs_Average_3DLCvar_JCScript
% function plot_Abs_Average_DLCvar_JCfun(VAR1, VAR2, VARname_str1, VARname_str2)
% Julien Catanese 9/23/2020

load('Average_ABS_DLC_LR.mat','VAR_str1','VAR_str2','VAR_str3',...
    'MeanABS1L','MeanABS1R','MeanABS2L','MeanABS2R','MeanABS3L','MeanABS3R',...
    'SemABS1L','SemABS1R','SemABS2L','SemABS2R','SemABS3L','SemABS3R',...
    'Thr1L','Thr1R','Thr2L','Thr2R','Thr3L','Thr3R');

%% LEFT vs RIGHT only plot 
figure, hold on, 

MeanABS1 = mean([MeanABS1L ],2);
MeanABS2 = mean([MeanABS2L ],2);
MeanABS3 = mean([MeanABS3L ],2);
SemABS1 = mean([SemABS1L ],2);
SemABS2 = mean([SemABS2L ],2);
SemABS3 = mean([SemABS3L ],2);
X=[1:1:100]';
plot(MeanABS1,'r','LineWidth',2);
plot(MeanABS2,'m','LineWidth',2);
plot(MeanABS3,'y','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
% plotshaded(X', [MeanABS1-(SemABS1*K) MeanABS1+(SemABS1*K)],'g');
% plotshaded(X', [MeanABS2-(SemABS2*K) MeanABS2+(SemABS2*K)],'c');
% plotshaded(X', [MeanABS3-(SemABS3*K) MeanABS3+(SemABS3*K)],'y');

MeanABS1 = mean([MeanABS1R ],2);
MeanABS2 = mean([MeanABS2R ],2);
MeanABS3 = mean([MeanABS3R ],2);
SemABS1 = mean([SemABS1R ],2);
SemABS2 = mean([SemABS2R ],2);
SemABS3 = mean([SemABS3R ],2);
X=[1:1:100]';
plot(MeanABS1,'b','LineWidth',2);
plot(MeanABS2,'c','LineWidth',2);
plot(MeanABS3,'g','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
% plotshaded(X', [MeanABS1-(SemABS1*K) MeanABS1+(SemABS1*K)],'g');
% plotshaded(X', [MeanABS2-(SemABS2*K) MeanABS2+(SemABS2*K)],'c');
% plotshaded(X', [MeanABS3-(SemABS3*K) MeanABS3+(SemABS3*K)],'y');
legend(VAR_str1 ,VAR_str2, VAR_str3  ,'GoCue');
title( ['Left Right DLC events' MouseID]);
% saveas(gcf, ['D:\DATA EMORY\JC_Analysis\Average_3DLCresults_LR' MouseID '_' VAR_str1 '_' VAR_str2 '_' VAR_str3 ],'png')
% saveas(gcf, ['D:\DATA EMORY\JC_Analysis\Average_3DLCresults_LR' MouseID '_' VAR_str1 '_' VAR_str2 '_' VAR_str3 ],'emf')
%% 

% MeanABS1 = mean([MeanABS1L MeanABS1R],2);
% MeanABS2 = mean([MeanABS2L MeanABS2R],2);
% MeanABS3 = mean([MeanABS3L MeanABS3R],2);
% 
% SemABS1 = mean([SemABS1L SemABS1R],2);
% SemABS2 = mean([SemABS2L SemABS2R],2);
% SemABS3 = mean([SemABS3L SemABS3R],2);
% 
% X=[1:1:100]';
% plot(MeanABS1,'g','LineWidth',2);
% plot(MeanABS2,'c','LineWidth',2);
% plot(MeanABS3,'y','LineWidth',2);
% xline(1.5*25,'k--','LineWidth',2);
% plotshaded(X', [MeanABS1-(SemABS1*K) MeanABS1+(SemABS1*K)],'g');
% plotshaded(X', [MeanABS2-(SemABS2*K) MeanABS2+(SemABS2*K)],'c');
% plotshaded(X', [MeanABS3-(SemABS3*K) MeanABS3+(SemABS3*K)],'y');
% 
% legend(VAR_str1 ,VAR_str2, VAR_str3  ,'GoCue');
% title( ['DLC events' MouseID]);
% 
% saveas(gcf, ['Average_3DLCresults_' MouseID '_' VAR_str1 '_' VAR_str2 '_' VAR_str3 ],'png')
% saveas(gcf, ['Average_3DLCresults_' MouseID '_' VAR_str1 '_' VAR_str2 '_' VAR_str3 ],'emf')
% saveas(gcf, ['D:\DATA EMORY\JC_Analysis\Average_3DLCresults_' MouseID '_' VAR_str1 '_' VAR_str2 '_' VAR_str3 ],'png')
% saveas(gcf, ['D:\DATA EMORY\JC_Analysis\Average_3DLCresults_' MouseID '_' VAR_str1 '_' VAR_str2 '_' VAR_str3 ],'emf')


