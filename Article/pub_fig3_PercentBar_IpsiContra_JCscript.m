% pub_fig_PercentBar_IpsiContra_JCscript
% plot a bar graph comparing percentages of cell type for contra and ipsi 
% by JC 2/19/2019 

close all 
load('listcell.mat')
figure, 

%% VMVL IPSI CONTRA
typecell = Tcombo.VMVL
ipsi_puff_cell = diff([TcLvR.Hcont_puf & typecell, TcLvR.Hipsi_puf & typecell ],1,2)>0; Nipsi_puff = sum(ipsi_puff_cell)
contra_puff_cell = diff([TcLvR.Hcont_puf & typecell , TcLvR.Hipsi_puf & typecell ],1,2)<0; Ncontra_puff = sum(contra_puff_cell)

ipsi_del_cell = diff([TcLvR.Hcont_del & typecell , TcLvR.Hipsi_del & typecell ],1,2)>0; Nipsi_del = sum(ipsi_del_cell)
contra_del_cell = diff([TcLvR.Hcont_del & typecell , TcLvR.Hipsi_del & typecell ],1,2)<0; Ncontra_del = sum(contra_del_cell)

ipsi_res_cell = diff([TcLvR.Hcont_res & typecell , TcLvR.Hipsi_res & typecell ],1,2)>0; Nipsi_res = sum(ipsi_res_cell)
contra_res_cell = diff([TcLvR.Hcont_res & typecell , TcLvR.Hipsi_res & typecell ],1,2)<0; Ncontra_res = sum(contra_res_cell)

ipsi_cell_tot = ipsi_puff_cell | ipsi_del_cell | ipsi_res_cell ; Nipsi_TOT = sum(ipsi_cell_tot)
contra_cell_tot = contra_puff_cell | contra_del_cell | contra_res_cell; Ncontra_TOT = sum(contra_cell_tot)

disp('type= VM/VL')
perc_ipsi = Nipsi_TOT/sum(typecell) ;
perc_contra = Ncontra_TOT/sum(typecell) ;
Ntot_typecell = sum(typecell) 
disp(['#ipsi=' num2str(Nipsi_TOT) ' (' num2str(round(perc_ipsi*100)) '%)'])
disp(['#contra=' num2str(Ncontra_TOT) ' (' num2str(round(perc_contra*100)) '%)'])

subplot(2,6,[1]), 
c = categorical({'ipsi','contra'});
bar(c,[perc_ipsi*100, perc_contra*100]), ylabel('%'), ylim([0 70]), 
title('VM')
DeltaContIps_VM = perc_contra - perc_ipsi; 

%% OPTO IPSI CONTRA
typecell =  (Tcombo.Opto_inib | Tcombo.Opto_exct) & Tcombo.VMVL ; 
ipsi_puff_cell = diff([TcLvR.Hcont_puf & typecell   , TcLvR.Hipsi_puf  & typecell ],1,2)>0; Nipsi_puff = sum(ipsi_puff_cell)
contra_puff_cell = diff([TcLvR.Hcont_puf  & typecell , TcLvR.Hipsi_puf  & typecell ],1,2)<0; Ncontra_puff = sum(contra_puff_cell)

ipsi_del_cell = diff([TcLvR.Hcont_del  & typecell  , TcLvR.Hipsi_del  & typecell  ],1,2)>0; Nipsi_del = sum(ipsi_del_cell)
contra_del_cell = diff([TcLvR.Hcont_del  & typecell , TcLvR.Hipsi_del  & typecell ],1,2)<0; Ncontra_del = sum(contra_del_cell)

ipsi_res_cell = diff([TcLvR.Hcont_res  & typecell , TcLvR.Hipsi_res  & typecell ],1,2)>0; Nipsi_res = sum(ipsi_res_cell)
contra_res_cell = diff([TcLvR.Hcont_res  & typecell , TcLvR.Hipsi_res  & typecell ],1,2)<0; Ncontra_res = sum(contra_res_cell)

ipsi_cell_tot = ipsi_puff_cell | ipsi_del_cell | ipsi_res_cell ; Nipsi_TOT = sum(ipsi_cell_tot)
contra_cell_tot = contra_puff_cell | contra_del_cell | contra_res_cell; Ncontra_TOT = sum(contra_cell_tot)

disp('type= OPTO')
perc_ipsi = Nipsi_TOT/sum(typecell) ;
perc_contra = Ncontra_TOT/sum(typecell) ;
Ntot_typecell = sum(typecell) 
disp(['#ipsi=' num2str(Nipsi_TOT) ' (' num2str(round(perc_ipsi*100)) '%)'])
disp(['#contra=' num2str(Ncontra_TOT) ' (' num2str(round(perc_contra*100)) '%)'])

subplot(2,6,[2]), 
c = categorical({'ipsi','contra'});
bar(c,[perc_ipsi*100, perc_contra*100]), ylabel('%'), ylim([0 70]), 
title('Opto')
DeltaContIps_opto = perc_contra - perc_ipsi; 

%% Puff IPSI CONTRA
disp('type= Puff')
typecell = Tcombo.PuffCell & Tcombo.VMVL; 
ipsi_puff_cell = diff([TcLvR.Hcont_puf  & typecell   , TcLvR.Hipsi_puf  & typecell ],1,2)>0; Nipsi_puff = sum(ipsi_puff_cell)
contra_puff_cell = diff([TcLvR.Hcont_puf  & typecell , TcLvR.Hipsi_puf  & typecell ],1,2)<0; Ncontra_puff = sum(contra_puff_cell)

ipsi_del_cell = diff([TcLvR.Hcont_del  & typecell  , TcLvR.Hipsi_del  & typecell  ],1,2)>0; Nipsi_del = sum(ipsi_del_cell)
contra_del_cell = diff([TcLvR.Hcont_del  & typecell , TcLvR.Hipsi_del  & typecell ],1,2)<0; Ncontra_del = sum(contra_del_cell)

ipsi_res_cell = diff([TcLvR.Hcont_res  & typecell , TcLvR.Hipsi_res  & typecell ],1,2)>0; Nipsi_res = sum(ipsi_res_cell)
contra_res_cell = diff([TcLvR.Hcont_res  & typecell , TcLvR.Hipsi_res  & typecell ],1,2)<0; Ncontra_res = sum(contra_res_cell)

ipsi_cell_tot = ipsi_puff_cell | ipsi_del_cell | ipsi_res_cell ; Nipsi_TOT = sum(ipsi_cell_tot)
contra_cell_tot= contra_puff_cell | contra_del_cell | contra_res_cell; Ncontra_TOT = sum(contra_cell_tot)

disp('type= Puff')
perc_ipsi = Nipsi_TOT/sum(typecell) ;
perc_contra = Ncontra_TOT/sum(typecell) ;
Ntot_typecell = sum(typecell) 
disp(['#ipsi=' num2str(Nipsi_TOT) ' (' num2str(round(perc_ipsi*100)) '%)'])
disp(['#contra=' num2str(Ncontra_TOT) ' (' num2str(round(perc_contra*100)) '%)'])

subplot(2,6,[3]), 
c = categorical({'ipsi','contra'});
bar(c,[perc_ipsi*100, perc_contra*100]), ylabel('%'), ylim([0 70]), 
title('Type1')
DeltaContIps_type1 = perc_contra - perc_ipsi; 

%% Delay IPSI CONTRA
disp('type= Delay')
typecell = Tcombo.Type2Cell & Tcombo.VMVL; 
ipsi_puff_cell = diff([TcLvR.Hcont_puf  & typecell   , TcLvR.Hipsi_puf  & typecell ],1,2)>0; Nipsi_puff = sum(ipsi_puff_cell)
contra_puff_cell = diff([TcLvR.Hcont_puf  & typecell , TcLvR.Hipsi_puf  & typecell ],1,2)<0; Ncontra_puff = sum(contra_puff_cell)

ipsi_del_cell = diff([TcLvR.Hcont_del  & typecell  , TcLvR.Hipsi_del  & typecell  ],1,2)>0; Nipsi_del = sum(ipsi_del_cell)
contra_del_cell = diff([TcLvR.Hcont_del  & typecell , TcLvR.Hipsi_del  & typecell ],1,2)<0; Ncontra_del = sum(contra_del_cell)

ipsi_res_cell = diff([TcLvR.Hcont_res  & typecell , TcLvR.Hipsi_res  & typecell ],1,2)>0; Nipsi_res = sum(ipsi_res_cell)
contra_res_cell = diff([TcLvR.Hcont_res  & typecell , TcLvR.Hipsi_res  & typecell ],1,2)<0; Ncontra_res = sum(contra_res_cell)

ipsi_cell_tot= ipsi_puff_cell | ipsi_del_cell | ipsi_res_cell ; Nipsi_TOT = sum(ipsi_cell_tot)
contra_cell_tot= contra_puff_cell | contra_del_cell | contra_res_cell; Ncontra_TOT = sum(contra_cell_tot)

disp('type= Delay')
perc_ipsi = Nipsi_TOT/sum(typecell) ;
perc_contra = Ncontra_TOT/sum(typecell) ;
Ntot_typecell = sum(typecell) 
disp(['#ipsi=' num2str(Nipsi_TOT) ' (' num2str(round(perc_ipsi*100)) '%)'])
disp(['#contra=' num2str(Ncontra_TOT) ' (' num2str(round(perc_contra*100)) '%)'])

subplot(2,6,[4]), 
c = categorical({'ipsi','contra'});
bar(c,[perc_ipsi*100, perc_contra*100]), ylabel('%'), ylim([0 70]), 
title('Type2')
DeltaContIps_type2 = perc_contra - perc_ipsi; 


%% Ramp IPSI CONTRA
disp('type= RAMP')
typecell = Tcombo.RespCell & Tcombo.VMVL; 

ipsi_puff_cell = diff([TcLvR.Hcont_puf  & typecell   , TcLvR.Hipsi_puf  & typecell ],1,2)>0; Nipsi_puff = sum(ipsi_puff_cell)
contra_puff_cell = diff([TcLvR.Hcont_puf  & typecell , TcLvR.Hipsi_puf  & typecell ],1,2)<0; Ncontra_puff = sum(contra_puff_cell)

ipsi_del_cell = diff([TcLvR.Hcont_del  & typecell  , TcLvR.Hipsi_del  & typecell  ],1,2)>0; Nipsi_del = sum(ipsi_del_cell)
contra_del_cell = diff([TcLvR.Hcont_del  & typecell , TcLvR.Hipsi_del  & typecell ],1,2)<0; Ncontra_del = sum(contra_del_cell)

ipsi_res_cell = diff([TcLvR.Hcont_res  & typecell , TcLvR.Hipsi_res  & typecell ],1,2)>0; Nipsi_res = sum(ipsi_res_cell)
contra_res_cell = diff([TcLvR.Hcont_res  & typecell , TcLvR.Hipsi_res  & typecell ],1,2)<0; Ncontra_res = sum(contra_res_cell)

ipsi_cell_tot = ipsi_puff_cell | ipsi_del_cell | ipsi_res_cell ; Nipsi_TOT = sum(ipsi_cell_tot)
contra_cell_tot= contra_puff_cell | contra_del_cell | contra_res_cell; Ncontra_TOT = sum(contra_cell_tot)

disp('type= RAMP')
perc_ipsi = Nipsi_TOT/sum(typecell) ;
perc_contra = Ncontra_TOT/sum(typecell) ;
Ntot_typecell = sum(typecell) 
disp(['#ipsi=' num2str(Nipsi_TOT) ' (' num2str(round(perc_ipsi*100)) '%)'])
disp(['#contra=' num2str(Ncontra_TOT) ' (' num2str(round(perc_contra*100)) '%)'])

subplot(2,6,[5]), 
c = categorical({'ipsi','contra'});
bar(c,[perc_ipsi*100, perc_contra*100]), ylabel('%'), ylim([0 70]), 
title('Type3')
DeltaContIps_type3 = perc_contra - perc_ipsi; 

%% Both IPSI CONTRA
disp('type= Both')
typecell = Tcombo.BothCell & Tcombo.VMVL; 
ipsi_puff_cell = diff([TcLvR.Hcont_puf  & typecell   , TcLvR.Hipsi_puf  & typecell ],1,2)>0; Nipsi_puff = sum(ipsi_puff_cell)
contra_puff_cell = diff([TcLvR.Hcont_puf  & typecell , TcLvR.Hipsi_puf  & typecell ],1,2)<0; Ncontra_puff = sum(contra_puff_cell)

ipsi_del_cell = diff([TcLvR.Hcont_del  & typecell  , TcLvR.Hipsi_del  & typecell  ],1,2)>0; Nipsi_del = sum(ipsi_del_cell)
contra_del_cell = diff([TcLvR.Hcont_del  & typecell , TcLvR.Hipsi_del  & typecell ],1,2)<0; Ncontra_del = sum(contra_del_cell)

ipsi_res_cell = diff([TcLvR.Hcont_res  & typecell , TcLvR.Hipsi_res  & typecell ],1,2)>0; Nipsi_res = sum(ipsi_res_cell)
contra_res_cell = diff([TcLvR.Hcont_res  & typecell , TcLvR.Hipsi_res  & typecell ],1,2)<0; Ncontra_res = sum(contra_res_cell)

ipsi_cell_tot= ipsi_puff_cell | ipsi_del_cell | ipsi_res_cell ; Nipsi_TOT = sum(ipsi_cell_tot)
contra_cell_tot= contra_puff_cell | contra_del_cell | contra_res_cell; Ncontra_TOT = sum(contra_cell_tot)

disp('type= Both')
perc_ipsi = Nipsi_TOT/sum(typecell) ;
perc_contra = Ncontra_TOT/sum(typecell) ;
Ntot_typecell = sum(typecell) 
disp(['#ipsi=' num2str(Nipsi_TOT) ' (' num2str(round(perc_ipsi*100)) '%)'])
disp(['#contra=' num2str(Ncontra_TOT) ' (' num2str(round(perc_contra*100)) '%)'])

subplot(2,6,[6]), 
c = categorical({'ipsi','contra'});
bar(c,[perc_ipsi*100, perc_contra*100]), ylabel('%'), 
ylim([0 70]), 
title('Type4')
DeltaContIps_type4 = perc_contra - perc_ipsi; 

%% Summary 
subplot(2,6,[7:12]), 
c = categorical({'Cell AllVM','Cell OptoVM','Type1(Puff)', 'Type2(Delay)', 'Type3(Resp)', 'Type4(Both)'});
bar(c,[DeltaContIps_VM*100, DeltaContIps_opto*100, DeltaContIps_type1*100, DeltaContIps_type2*100, DeltaContIps_type3*100, DeltaContIps_type4*100]), 
ylabel('%Contra - %Ipsi'), 
ylim([-30 +15]), 
title('Diff Contra-Ipsi')



