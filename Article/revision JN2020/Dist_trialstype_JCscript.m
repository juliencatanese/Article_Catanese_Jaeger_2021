
clear all, close all,
%% To run a loop over folder
AnalysisFolder = 'D:\DATA EMORY\JC_Analysis\';
MouseID = 'JCVGAT';
FolderID = dir([AnalysisFolder MouseID '*\*task*']);

DistPerc_cor = [];    DistPerc_imp = []; DistPerc_omi = []; NTOTAL = 0;

for nf=1:max(size(FolderID));
    cd([FolderID(nf).folder '\' FolderID(nf).name]);
    disp([num2str(nf) '/' num2str(max(size(FolderID))) '   MouseID =' MouseID  '   FolderID = ' FolderID(nf).name]);
    
    load('Ntrial_type.mat');
    ntot = trial.Ntrial;
    idx_cor = [trial.idx_correct_R  trial.idx_correct_L];
    idx_imp = [trial.idx_errorDelay_PL_CL trial.idx_errorDelay_PR_CL trial.idx_errorDelay_PL_CR trial.idx_errorDelay_PR_CR];
    idx_omi = [trial.idx_NoLick];
    
    DistPerc_cor = [DistPerc_cor (idx_cor/ntot)*100];
    DistPerc_imp = [DistPerc_imp (idx_imp/ntot)*100];
    DistPerc_omi = [DistPerc_omi (idx_omi/ntot)*100];
    NTOTAL = NTOTAL + ntot;
end
%%
cd(AnalysisFolder)
nbin  = 10
NormT = {'count','probability','cumcount'};
for ii= 1:3
    figure, hold on;
    histogram(DistPerc_cor, nbin, 'Normalization', NormT{ii});
    histogram(DistPerc_imp, nbin, 'Normalization', NormT{ii});
    histogram(DistPerc_omi, nbin, 'Normalization', NormT{ii});
    legend('cor','imp','omi','Location','best')
    xlabel('percentage of total session lenght')
    ylabel(NormT{ii})
    title([ '#' num2str(nbin) 'bins  ' '#' num2str(NTOTAL) 'trials  '  '#' num2str(nf) 'Sessions  '  '#5mice'])
    saveas(gcf, ['ArticleFig_DIST_TRIALTYPE_' NormT{ii} '_'  num2str(nbin) '#bins' ],'png');
    saveas(gcf, ['ArticleFig_DIST_TRIALTYPE_' NormT{ii} '_'  num2str(nbin) '#bins' ],'emf');
end

%%
NormT = {'PercTrialType'};  ii=1;
A= sort(DistPerc_cor);
B= sort(DistPerc_imp);
C= sort(DistPerc_omi);
BINlist= [3 4 5 10]
for nn=1:max(size(BINlist))
    nbin= BINlist(nn)
    Wbin= (100/nbin)
    
    AA=[]; BB=[]; CC=[];AAA=[]; BBB=[];CCC=[];
    for jj=1:nbin
        percValUp= (100/nbin)*jj
        percValDown= (100/nbin)*(jj-1)
        LABEL{jj} = num2str(percValUp)
        
        AA(jj) = sum(A<percValUp & A>=percValDown)
        BB(jj) = sum(B<percValUp & B>=percValDown)
        CC(jj) = sum(C<percValUp & C>=percValDown)
        
        
        Nbintot(jj) = AA(jj) + BB(jj) + CC(jj)
        AAA(jj)= 100*(AA(jj)/Nbintot(jj))
        BBB(jj)= 100*(BB(jj)/Nbintot(jj))
        CCC(jj)= 100*(CC(jj)/Nbintot(jj))
        
    end
    
    f1 = figure; ax1 = axes('Parent',f1); hold(ax1,'on');
    bar([AAA'  BBB' CCC'],'stacked')
    % set(ax1,'XTick',[1:1:nbin],'XTickLabel',{'20','40','60','80','100'});
    set(ax1,'XTick',[1:1:nbin],'XTickLabel',LABEL);
    xlabel('percentage of total session lenght (%)')
    ylabel('percentage of trials types (%)')
    legend('cor','imp','omi','Location','best')
    title([ '#' num2str(nbin) 'bins  ' '#' num2str(NTOTAL) 'trials  '  '#' num2str(nf) 'Sessions  '  '#5mice'])
    saveas(gcf, ['ArticleFig_DIST_TRIALTYPE_' NormT{ii} '_'  num2str(nbin) '#bins' ],'png');
    saveas(gcf, ['ArticleFig_DIST_TRIALTYPE_' NormT{ii} '_'  num2str(nbin) '#bins' ],'emf');
end