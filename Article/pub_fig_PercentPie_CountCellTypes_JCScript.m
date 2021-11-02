% pub_fig_PercentPie_CountCellTypes_JCScript
% save eps figure of the percentage of cell types (Pie)
% by JC 2/19/2019

close all
clear all
load ('listcell.mat')

%% Ephys Type in VMVL cells
Npuff = sum(Tcombo.PuffCell & Tcombo.VMVL); disp(['#VMVL Puff Cell = ' num2str(Npuff)])
Nboth= sum(Tcombo.BothCell & Tcombo.VMVL); disp(['#VMVL Both Cell = ' num2str(Nboth)])
Nresp= sum(Tcombo.RespCell & Tcombo.VMVL); disp(['#VMVL Resp Cell = ' num2str(Nresp)])

Ntype1 = sum(Tcombo.Type1Cell & Tcombo.VMVL); disp(['#VMVL Type1 Cell = ' num2str(Ntype1)])
Ntype2 = sum(Tcombo.Type2Cell & Tcombo.VMVL); disp(['#VMVL Type2 Cell = ' num2str(Ntype2)])
Ntype3 = sum(Tcombo.Type3Cell & Tcombo.VMVL); disp(['#VMVL Type3 Cell = ' num2str(Ntype3)])
other = sum(Tcombo.NoSigCell  & Tcombo.VMVL); disp(['#VMVL NoSIg Cell = ' num2str(other)])

Ntot=sum(Tcombo.VMVL);

labels = {['Type1 (Sample ' num2str(round(Npuff/Ntot*100)) '%)  '],...
    ['Type2 (Delay ' num2str(round(Ntype2/Ntot*100)) '%) '],...
    ['Type3 (Resp ' num2str(round(Nresp/Ntot*100)) '%) ' ],...
    [ 'Type4 (Both ' num2str(round(Nboth/Ntot*100)) '%)' ],...
    [ 'Type5 (Other ' num2str(round(other/Ntot*100)) '%) ' ]};
explode = [0 0 1 1 0];
figure,
pie([Npuff/Ntot, Ntype2/Ntot, Nresp/Ntot, Nboth/Ntot, other/Ntot], explode, labels)
title('Ephys CellTypes VMVL NONopto')
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie1_EphysVM'],'eps')
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie1_EphysVM'],'emf')
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie1_EphysVM'],'jpg')


%% Opto + VMVL + Type
% OPTO all
Ntot_opto= sum(Tcombo.Opto_post_sess & Tcombo.VMVL);
Nopto_inhib= sum(Tcombo.Opto_inib & Tcombo.VMVL);
Nopto_excit= sum(Tcombo.Opto_exct & Tcombo.VMVL);
Nopto_nonresp = Ntot_opto - Nopto_inhib - Nopto_excit;

labels = {['opto inhib (' num2str(round(Nopto_inhib/Ntot_opto*100)) '%)  '],...
    ['opto excite (' num2str(round(Nopto_excit/Ntot_opto*100)) '%) '],...
    ['opto nonResp (' num2str(round(Nopto_nonresp/Ntot_opto*100)) '%) ' ]};

figure,
explode=[0 0 1]
pie([Nopto_inhib/Ntot_opto, Nopto_excit/Ntot_opto, Nopto_nonresp/Ntot_opto ], explode, labels)
title('OPTO SNr-VMVL')
disp(['Saving : D:\JC_Figures\Fig_Article\Pie2_optoEphysVM'])
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie2_optoEphysVM'],'eps')
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie2_optoEphysVM'],'emf')
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie2_optoEphysVM'],'jpg')

%% OPTO INHIB
Nopto_inhib

PUFF_VMVL_Opto_Inib = sum(Tcombo.VMVL & Tcombo.Opto_inib & Tcombo.PuffCell) ;  disp(['# PUFF_VMVL_Opto_Inib = ' num2str(PUFF_VMVL_Opto_Inib)])
BOTH_VMVL_Opto_Inib = sum(Tcombo.VMVL & Tcombo.Opto_inib & Tcombo.BothCell) ;  disp(['# BOTH_VMVL_Opto_Inib = ' num2str(BOTH_VMVL_Opto_Inib)])
RESP_VMVL_Opto_Inib = sum(Tcombo.VMVL & Tcombo.Opto_inib & Tcombo.RespCell) ; disp(['# RESP_VMVL_Opto_Inib = ' num2str(RESP_VMVL_Opto_Inib)])
TYPE1_VMVL_Opto_Inib = sum(Tcombo.VMVL & Tcombo.Opto_inib & Tcombo.Type1Cell);  disp(['# TYPE1_VMVL_Opto_Inib = ' num2str(TYPE1_VMVL_Opto_Inib)])
TYPE2_VMVL_Opto_Inib = sum(Tcombo.VMVL & Tcombo.Opto_inib & Tcombo.Type2Cell);  disp(['# TYPE2_VMVL_Opto_Inib = ' num2str(TYPE2_VMVL_Opto_Inib)])
TYPE3_VMVL_Opto_Inib = sum(Tcombo.VMVL & Tcombo.Opto_inib & Tcombo.Type3Cell);  disp(['# TYPE3_VMVL_Opto_Inib = ' num2str(TYPE3_VMVL_Opto_Inib)])
NOSIG_VMVL_Opto_Inib = sum(Tcombo.VMVL & Tcombo.Opto_inib & Tcombo.NoSigCell);  disp(['# NOSIG_VMVL_Opto_Inib = ' num2str(NOSIG_VMVL_Opto_Inib)])

labels = {['Type1 (Sample ' num2str(round(PUFF_VMVL_Opto_Inib/Nopto_inhib*100)) '%)  '],...
    ['Type2 (Delay ' num2str(round(TYPE2_VMVL_Opto_Inib/Nopto_inhib*100)) '%) '],...
    ['Type3 (Resp ' num2str(round(RESP_VMVL_Opto_Inib/Nopto_inhib*100)) '%) ' ],...
    [ 'Type4 (Both ' num2str(round(BOTH_VMVL_Opto_Inib/Nopto_inhib*100)) '%)' ],...
    [ 'Type5 (Other ' num2str(round(NOSIG_VMVL_Opto_Inib/Nopto_inhib*100)) '%) ' ]};
explode = [0 0 1 1 0];
figure,
pie([PUFF_VMVL_Opto_Inib/Nopto_inhib, TYPE2_VMVL_Opto_Inib/Nopto_inhib, RESP_VMVL_Opto_Inib/Nopto_inhib, BOTH_VMVL_Opto_Inib/Nopto_inhib, NOSIG_VMVL_Opto_Inib/Nopto_inhib], explode, labels)
title(['OPTO Inhib SNr-VMVL Ephys CellTypes NtotInhib=' num2str(Nopto_inhib) ])
disp(['Saving : D:\JC_Figures\Fig_Article\Pie3_optoINHIB'])
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie3_optoINHIB'],'eps')
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie3_optoINHIB'],'emf')
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie3_optoINHIB'],'jpg')
%
%% OPTO EXCIT
Nopto_excit

PUFF_VMVL_Opto_Excit = sum(Tcombo.VMVL & Tcombo.Opto_exct & Tcombo.PuffCell) ;  disp(['# PUFF_VMVL_Opto_Excit = ' num2str(PUFF_VMVL_Opto_Excit)])
BOTH_VMVL_Opto_Excit = sum(Tcombo.VMVL & Tcombo.Opto_exct & Tcombo.BothCell) ;  disp(['# BOTH_VMVL_Opto_Excit  = ' num2str(BOTH_VMVL_Opto_Excit )])
RESP_VMVL_Opto_Excit = sum(Tcombo.VMVL & Tcombo.Opto_exct & Tcombo.RespCell) ;  disp(['# RESP_VMVL_Opto_Excit = ' num2str(RESP_VMVL_Opto_Excit)])
TYPE1_VMVL_Opto_Excit = sum(Tcombo.VMVL & Tcombo.Opto_exct & Tcombo.Type1Cell) ;  disp(['# TYPE1_VMVL_Opto_Excit = ' num2str(TYPE1_VMVL_Opto_Excit)])
TYPE2_VMVL_Opto_Excit = sum(Tcombo.VMVL & Tcombo.Opto_exct & Tcombo.Type2Cell) ;  disp(['# TYPE2_VMVL_Opto_Excit = ' num2str(TYPE2_VMVL_Opto_Excit)])
TYPE3_VMVL_Opto_Excit = sum(Tcombo.VMVL & Tcombo.Opto_exct & Tcombo.Type3Cell) ;  disp(['# TYPE3_VMVL_Opto_Excit = ' num2str(TYPE3_VMVL_Opto_Excit)])
NOSIG_VMVL_Opto_Excit = sum(Tcombo.VMVL & Tcombo.Opto_exct & Tcombo.NoSigCell) ;  disp(['# NOSIG_VMVL_Opto_Excit = ' num2str(NOSIG_VMVL_Opto_Excit)])

labels = {['Type1 (Sample ' num2str(round(PUFF_VMVL_Opto_Excit/Nopto_excit*100)) '%)  '],...
    ['Type2 (Delay ' num2str(round(TYPE2_VMVL_Opto_Excit/Nopto_excit*100)) '%) '],...
    ['Type3 (Resp ' num2str(round(RESP_VMVL_Opto_Excit/Nopto_excit*100)) '%) ' ],...
    [ 'Type4 (Both ' num2str(round(BOTH_VMVL_Opto_Excit/Nopto_excit*100)) '%)' ],...
    [ 'Type5 (Other ' num2str(round(NOSIG_VMVL_Opto_Excit/Nopto_excit*100)) '%) ' ]};
explode = [0 0 1 1 0];
figure,
pie([Npuff/Ntot, Ntype2/Ntot, Nresp/Ntot, Nboth/Ntot, other/Ntot], explode, labels)
title(['OPTO excit SNr-VMVL Ephys CellTypes NtotExcit=' num2str(Nopto_excit) ])
disp(['Saving : D:\JC_Figures\Fig_Article\Pie4_optoEXCIT'])
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie4_optoEXCIT'],'eps')
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie4_optoEXCIT'],'emf')
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie4_optoEXCIT'],'jpg')

