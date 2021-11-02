% pub_fig2_OptoPie_JCscript
% plot a pie-chart %cell modulated by POST OPTO STIM (+/-/0)
% Written by Julien Catanese 2/19/2019
% Updated JC 3/27/2019: only keep the Opto part. 

% Count OPTO responsive cell in VMVL
Ntot_opto= sum(Tcombo.Opto_post_sess & Tcombo.VMVL);
Nopto_inhib= sum(Tcombo.Opto_inib & Tcombo.VMVL);
Nopto_excit= sum(Tcombo.Opto_exct & Tcombo.VMVL);
Nopto_nonresp = Ntot_opto - Nopto_inhib - Nopto_excit;

% labels for pie legend OPTO responsive cell in VMVL
labels = {['opto inhib (' num2str(round(Nopto_inhib/Ntot_opto*100)) '%)  '],...
    ['opto excite (' num2str(round(Nopto_excit/Ntot_opto*100)) '%) '],...
    ['opto nonResp (' num2str(round(Nopto_nonresp/Ntot_opto*100)) '%) ' ]};

% Plot OPTO responsive cell in VMVL
figure,
explode=[0 0 1]
pie([Nopto_inhib/Ntot_opto, Nopto_excit/Ntot_opto, Nopto_nonresp/Ntot_opto ], explode, labels)
title('OPTO SNr-VMVL')
disp(['Saving : D:\JC_Figures\Fig_Article\Pie2_optoEphysVM'])
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie2_optoEphysVM'],'fig')
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie2_optoEphysVM'],'emf')
saveas(gcf, ['D:\JC_Figures\Fig_Article\Pie2_optoEphysVM'],'jpg')
