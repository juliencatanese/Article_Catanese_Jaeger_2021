% pub_tableFinal_Tcombo_JCscript
% Combine all the main info from all previous tables 
% by JC 2/19/2019

load listcell.mat
Topto.H_opto_post(isnan(Topto.H_opto_post))=0;
Ncell_opto =  sum(Topto.H_opto_post)
Ncell_opto_inib = sum(Topto.H_opto_post & (Topto.MeanNspk_optoON_post<Topto.MeanNspk_optoOFF_post))
Ncell_opto_exct = sum(Topto.H_opto_post & (Topto.MeanNspk_optoON_post>Topto.MeanNspk_optoOFF_post))

VMVL      = Tcoord.VM + Tcoord.VL;

AreaID = Tcoord.AreaID;  

Opto_inib = Topto.H_opto_post & (Topto.MeanNspk_optoON_post<Topto.MeanNspk_optoOFF_post);
Opto_exct = Topto.H_opto_post & (Topto.MeanNspk_optoON_post>Topto.MeanNspk_optoOFF_post);
Opto_post_sess = Topto.OPTO_POST;
Opto_post_Nstim = Topto.OPTO_POST_Nstim;

PuffCell  = (Tephys.H_puf & Tephys.H_del & ~Tephys.H_res) + (Tephys.H_puf & ~Tephys.H_del & ~Tephys.H_res)  ;
BothCell  = (Tephys.H_puf & ~Tephys.H_del & Tephys.H_res) + (Tephys.H_puf & Tephys.H_del & Tephys.H_res)  ;
RespCell  = (~Tephys.H_puf & Tephys.H_del & Tephys.H_res) + (~Tephys.H_puf & ~Tephys.H_del & Tephys.H_res) ;

Type1Cell  = (Tephys.H_puf  & ~Tephys.H_del  & ~Tephys.H_res) ;
Type2Cell  = (~Tephys.H_puf & Tephys.H_del  & ~Tephys.H_res);
Type3Cell  = (~Tephys.H_puf & ~Tephys.H_del & Tephys.H_res)  ;
NoSigCell = (~Tephys.H_puf & ~Tephys.H_del & ~Tephys.H_res)  ;

ipsi_puff_cell = diff([TcLvR.Hcont_puf & VMVL , TcLvR.Hipsi_puf & VMVL ],1,2)>0; Nipsi_puff = sum(ipsi_puff_cell);
contra_puff_cell = diff([TcLvR.Hcont_puf & VMVL , TcLvR.Hipsi_puf & VMVL ],1,2)<0; Ncontra_puff = sum(contra_puff_cell);

ipsi_del_cell = diff([TcLvR.Hcont_del & VMVL , TcLvR.Hipsi_del & VMVL ],1,2)>0; Nipsi_del = sum(ipsi_del_cell);
contra_del_cell = diff([TcLvR.Hcont_del & VMVL , TcLvR.Hipsi_del & VMVL ],1,2)<0; Ncontra_del = sum(contra_del_cell);

ipsi_res_cell = diff([TcLvR.Hcont_res & VMVL , TcLvR.Hipsi_res & VMVL ],1,2)>0; Nipsi_res = sum(ipsi_res_cell);
contra_res_cell = diff([TcLvR.Hcont_res & VMVL , TcLvR.Hipsi_res & VMVL ],1,2)<0; Ncontra_res = sum(contra_res_cell);

ipsi_cell = ipsi_puff_cell | ipsi_del_cell | ipsi_res_cell ; Nipsi_TOT = sum(ipsi_cell);
contra_cell = contra_puff_cell | contra_del_cell | contra_res_cell; Ncontra_TOT = sum(contra_cell);

Tcombo= [];
Tcombo= addvars(Tcoord(:,1:6), VMVL, AreaID, Opto_inib, Opto_exct, ipsi_cell, contra_cell, RespCell, PuffCell, BothCell, Type1Cell, Type2Cell, Type3Cell, NoSigCell, Opto_post_sess, Opto_post_Nstim, listcell(:,2:3));
Tcombo(1,:);
save('D:\JC_Analysis\listcell.mat','Tcombo', '-append');
disp('Tcombo saved')

% Rename listcell
if parfig.center_evt=='Delay'
    [Sdel]=copyfile('D:\JC_Analysis\listcell.mat', 'D:\JC_Analysis\listcell_DelayCenter.mat')
elseif parfig.center_evt=='GoCue'
    [Scor]=copyfile('D:\JC_Analysis\listcell.mat', 'D:\JC_Analysis\listcell_GoCueCenter.mat')
end

