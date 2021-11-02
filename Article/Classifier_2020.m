% pub_table_MasterScript_JCscript
% by JC 5/6/2020

%% Fig7  Classifier
% DELAY epoch, GOCue Centered, cor vs omi, cor vs imp
% Computing the Classifier10FOld _Mall.mat (saved in JC_Analysis)
close all;
mypath = 'D:\DATA EMORY\JC_Analysis'
cd(mypath)
parfig.pathinit = mypath
SaveFigFolder = 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig7_classifier'
clearvars -except mypath parfig SaveFigFolder
parfig.ControlShuffle=0;
rng('shuffle');
load('listcell3.mat'); 
load('Tfig1_VMopto.mat');

parfig.Nrepeat = 10;
parfig.typeClass='FoldXVal' ; %     typeClass='%HoldOut'
parfig.learner= 'logistic';%, 'svm'}
parfig.Nfold = 10;

parfig.center_evt = 'GoCue';% center ('Delay' ; 'GoCue' ;'APuff'; 'Licks')
parfig.epoch = 'delay'

if parfig.epoch  =='delay'
    parfig.pre = 750 ; % define how much time before zero (in ms)
    parfig.post= 0;
elseif parfig.epoch == 'Apuff'
    parfig.pre = 1500 ; % define how much time before zero (in ms)
    parfig.post= -750;
end

%% Fig7A All raw traces per sessions (very long to recompute)
% VMVL all
celltype = 'VMVL'
idxcell = Tfig1_VMopto.VMVL
parfig.minNbcell = 20;
parfig.ControlShuffle=0;

trial_type_1 = 'cor'
trial_type_2 = 'imp'
pub_fig7_Classifier_JCscript % very long

%% Average for distinct population
% PRE-LICK epoch, Lick Centered, cor vs imp
% PLOTING Mean +/- 1SEM (shaded)
SaveFigFolder = 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig7_classifier'
clearvars -except mypath parfig Tfig1_VMopto listcell SaveFigFolder jj
figure(10) ; close (10); figure(10)
parfig.center_evt = 'Licks';    % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.epoch = 'delay'

for ii=2:3
    if ii==1
        celltype = 'NOL'
        parfig.minNbcell = 15;
        col = 'k'
    elseif ii==2
        celltype = 'SvTh-'
        parfig.minNbcell = 9;
        col = 'c'
    elseif ii==3
        celltype = 'ShuffleControl'
        parfig.minNbcell = 20;
        col = 'r'
    end
    
    trial_type_1 = 'cor'
    trial_type_2 = 'imp';
    
    load(['Class10Fold_Mall_' trial_type_1 'v' trial_type_2  '_' parfig.epoch '_' parfig.center_evt '_' celltype '.mat'], 'Mall')
    figure (10), hold on,
    K = 1;

    Mm=mean(Mall)
    NSessFinal = size(Mall,1)
    Mstd = std(Mall)/sqrt(NSessFinal)
    XX=[1:1:parfig.minNbcell];
    plot( XX , Mm, col, 'LineWidth', 3)
    plotshaded(XX, [Mm-K*Mstd ; Mm+K*Mstd], col)
end

hold on,
plot( XX , ones(1,parfig.minNbcell)/2, 'k--')
ylim([0.2 0.8]);
title(['Mean Class10Fold ' trial_type_1 ' vs ' trial_type_2  ' ' parfig.epoch  ' ' parfig.center_evt ' '  celltype ])
legend('NOL', '', 'SvTh-', '','Shuffle','')

