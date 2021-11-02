% Select_Sess_trialTrsh_JCscript

%%
AnalysisFolder = 'C:\Users\catan\Documents\EMORY\JC_Analysis\';
MouseID = 'JCVGAT';
FolderID = dir([AnalysisFolder MouseID '*\*task*']);

bool_select = []; 
List_Sess = [];
for nf=1:max(size(FolderID))
    cd([FolderID(nf).folder '\' FolderID(nf).name])
    disp([num2str(nf) '/' num2str(max(size(FolderID))) '   MouseID =' MouseID  '   FolderID = ' FolderID(nf).name])
    
    List_Sess = [List_Sess ; FolderID(nf).name(1:12)]
    
    Count_tr_types_JC_Script;
    ALL_opto_tr = trial.Nb_all_opto
    imp_opto_tr = trial.Nb_errorDelay_PL_CL_opto + trial.Nb_errorDelay_PL_CR_opto +  trial.Nb_errorDelay_PR_CL_opto + trial.Nb_errorDelay_PR_CR_opto
    omi_opto_tr = trial.Nb_NoLick_opto
    
    if imp_opto_tr >= 3 & omi_opto_tr >= 5
        bool_select = [bool_select 1]
    else
        bool_select = [bool_select 0]
    end
    
end

List_Sess_all = List_Sess
idx_Sess = find(bool_select)
List_Sess(idx_Sess,:)

for ii=1:length(idx_Sess)
boola(:,ii)= Tfig2_cor.nSess==idx_Sess(ii);
end
bool_sel = sum(boola,2);

cd (AnalysisFolder) 