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
idxL = trial.idx_correct_L; ntrL=size(idxL,2);
idxR = trial.idx_correct_R; ntrR=size(idxR,2);

Thr1L = mean(mean(VAR1(10:30,idxL)));Thr1R = mean(mean(VAR1(10:30,idxR)));
Thr2L = mean(mean(VAR2(10:30,idxL)));Thr2R = mean(mean(VAR2(10:30,idxR)));
Thr3L = mean(mean(VAR3(10:30,idxL)));Thr3R = mean(mean(VAR3(10:30,idxR)));
MeanABS1L = mean(abs(-1*(VAR1(:,idxL)-Thr1L)),2); MeanABS1R = mean(abs(VAR1(:,idxR)-Thr1R),2);
MeanABS2L = mean(abs(-1*(VAR2(:,idxL)-Thr2L)),2); MeanABS2R = mean(abs(VAR2(:,idxR)-Thr2R),2);
MeanABS3L = mean(abs(-1*(VAR3(:,idxL)-Thr3L)),2); MeanABS3R = mean(abs(VAR3(:,idxR)-Thr3R),2);

SemABS1L = (std(abs(-1*(VAR1(:,idxL)-Thr1L))',1)/sqrt(ntrL))'; SemABS1R = (std(abs(VAR1(:,idxR)-Thr1R)',1)/sqrt(ntrR))';
SemABS2L = (std(abs(-1*(VAR2(:,idxL)-Thr2L))',1)/sqrt(ntrL))'; SemABS2R = (std(abs(VAR2(:,idxR)-Thr2R)',1)/sqrt(ntrR))';
SemABS3L = (std(abs(-1*(VAR3(:,idxL)-Thr3L))',1)/sqrt(ntrL))'; SemABS3R = (std(abs(VAR3(:,idxR)-Thr3R)',1)/sqrt(ntrR))';

save('Average_ABS_DLC_LR.mat','VAR_str1','VAR_str2','VAR_str3',...
    'MeanABS1L','MeanABS1R','MeanABS2L','MeanABS2R','MeanABS3L','MeanABS3R',...
    'SemABS1L','SemABS1R','SemABS2L','SemABS2R','SemABS3L','SemABS3R',...
    'Thr1L','Thr1R','Thr2L','Thr2R','Thr3L','Thr3R');

% MeanABS1 = mean([MeanABS1L MeanABS1R],2);
% MeanABS2 = mean([MeanABS2L MeanABS2R],2);
% MeanABS3 = mean([MeanABS3L MeanABS3R],2);
% 
% SemABS1 = mean([SemABS1L SemABS1R],2);
% SemABS2 = mean([SemABS2L SemABS2R],2);
% SemABS3 = mean([SemABS3L SemABS3R],2);

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