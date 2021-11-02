
%% INSIDE pub_MAsterScript_JCscript=
%Fig7  plot Classifier Left vs Right during Delay
clear all,
close all,
pause(0.1);
mypath = 'C:\Users\catan\Documents\EMORY\JC_Analysis'
cd(mypath)

FolderID = dir([mypath '\JCVGAT*\*taskopto*'])
All_RTcor = [];
All_RTocC = [];
for nf=1: max(size(FolderID))
    FileLocation = [FolderID(nf).folder '\' FolderID(nf).name]
    MouseID= FolderID(nf).name(1:6)
    Day= FolderID(nf).name(8:12)
    cd(FileLocation)
    
    for jj=1:4
        if jj==1
            parfig.center_evt = 'GoCue';  Evt1= 'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
            parfig.trial_type = {'cor'}; trialtype1= parfig.trial_type{1}; tr1= 'cor';
        elseif jj==2
            parfig.center_evt = 'GoCue';  Evt2= 'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
            parfig.trial_type = {'ocC'}; trialtype1= parfig.trial_type{1}; tr2= 'ocC';
        elseif jj==3
            parfig.center_evt = 'Licks';  Evt3= 'Licks';  % center for SDF ('Delay' ; 'GOcue' ; Licks)
            parfig.trial_type = {'cor'}; trialtype1= parfig.trial_type{1}; tr3= 'cor';
        elseif jj==4
            parfig.center_evt = 'Licks';  Evt4= 'Licks'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
            parfig.trial_type = {'ocC'}; trialtype1= parfig.trial_type{1}; tr4= 'ocC';
        end
        
        % GET trigtimes : times of the psth_trig_evt is sec   (1 vector of times)
        [trigtimes] = Get_trigtimes(parfig.center_evt, parfig.trial_type, FileLocation); AA=trigtimes;  % trigtimes in sec
        
        if jj==1
            T1 = trigtimes;
        elseif jj==2
            T2 = trigtimes;
        elseif jj==3
            T3 = trigtimes;
        elseif jj==4
            T4 = trigtimes;
        end
    end
    
    %% display results
    % Verify the Numbers
    load([FileLocation '/Ntrial_type.mat'])
    trial
    Nb_cor_optoOFF =  trial.Nb_correct_L + trial.Nb_correct_R
    Nb_cor_optoON  =  trial.Nb_correct_L_opto + trial.Nb_correct_R_opto
    
    Evt1
    tr1
    Nbtimes_T1 = max(size(T1))
    
    Evt2
    tr2
    Nbtimes_T2 = max(size(T2))
    
    Evt3
    tr3
    Nbtimes_T3 = max(size(T3))
    
    Evt4
    tr4
    Nbtimes_T4 = max(size(T4))
    
    %% Reaction Time for each session
    RT_cor = T3-T1;
    RT_ocC = T4-T2;
    
    RT_cor(find(RT_cor<0))=[];
    RT_ocC(find(RT_ocC<0))=[];
    
    av_RT_cor = mean(RT_cor)
    av_RT_ocC = mean(RT_ocC)
    
    std_RT_cor = std(RT_cor)
    std_RT_ocC = std(RT_ocC)
    
    %% Global reaction time Distribution
    % Make one vector containing RT of each sessions
    All_RTcor = [All_RTcor RT_cor];
    All_RTocC = [All_RTocC RT_ocC];
    
    %%
    figure (nf), pause(0.01);
    hold on, histogram(RT_cor, 25), pause(0.01);
    hold on, histogram(RT_ocC, 25), pause(0.01);
    legend ('cor optoOFF','cor optoON')
    title ([MouseID Day])
    xlim ([0 2])
    pause(0.01)
    
end
%% PLOT OVERALL DISTRIBUTION
figure (100), close 100
figure (100), pause(0.01);
xlim ([0 2])
hold on, histogram(All_RTcor, 50), pause(0.01);
hold on, histogram(All_RTocC, 25), pause(0.01);
legend ('cor optoOFF','cor optoON')
title ([ '#' num2str(max(size([All_RTcor All_RTocC]))) ' trials' ' 11 sessions 4 mice ' ])
xlim ([0 2])
%%
figure,
boxplot([All_RTcor(1,222:221+222)' All_RTocC'], 'labels', {'cor optoOFF','cor optoON'} )
%     ranksum test for equal medians:
%     The Wilcoxon rank sum test is equivalent to the Mann-Whitney U test.
%     P = ranksum(X,Y) performs a two-sided rank sum test of the hypothesis that two independent samples,
%     in the vectors X and Y, come from distributions with equal medians, and returns the p-value from the test.
X= All_RTcor;
Y= All_RTocC;
[P,H] = ranksum(X,Y,'tail','both')
title(['Wilcoxon test P=' num2str(P)])
pause(0.01)
%% SMALL VALUE
close all
RTvalue = 0.25 %sec
smallRTcor= All_RTcor(find(All_RTcor<RTvalue));
smallRTocC= All_RTocC(find(All_RTocC<RTvalue));

% DIST
figure (200), close 200
figure (200), pause(0.01);
xlim ([0 RTvalue])
hold on, histogram(smallRTcor, 50), pause(0.01);
hold on, histogram(smallRTocC, 25), pause(0.01);
legend ('cor optoOFF','cor optoON')
title ([ '#' num2str(max(size([smallRTcor smallRTocC]))) ' trials with RT<' num2str(RTvalue) 'sec'])

% Average STAT
figure,
boxplot([smallRTcor(1:max(size(smallRTocC)))' smallRTocC' ], 'labels', {'cor optoOFF','cor optoON'} )
%     ranksum test for equal medians:
%     The Wilcoxon rank sum test is equivalent to the Mann-Whitney U test.
%     P = ranksum(X,Y) performs a two-sided rank sum test of the hypothesis that two independent samples,
%     in the vectors X and Y, come from distributions with equal medians, and returns the p-value from the test.
X= smallRTcor;
Y= smallRTocC;
[P,H] = ranksum(X,Y,'tail','left')
title(['Wilcoxon test P=' num2str(P)])
pause(0.01)
