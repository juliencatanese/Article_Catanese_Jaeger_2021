% pub_table6_Tcombo_JCscript
% Combine all the main info from all previous tables 
% Save additional Table in listcell.mat: Tcombo, Tcombo_z, Tcombo_zipsi...  
% written by JC 02/19/2019
% last updated JC 03/31/19: "add Opto_post_sess to all Tz"

load listcell.mat
%% From T.coord 
VMVL      = Tcoord.VM + Tcoord.VL;
AreaID = Tcoord.AreaID;  

%% From T.opto 
Topto.H_opto_post(isnan(Topto.H_opto_post))=0;
Opto_inib = Topto.H_opto_post & (Topto.MeanNspk_optoON_post<Topto.MeanNspk_optoOFF_post);
Opto_exct = Topto.H_opto_post & (Topto.MeanNspk_optoON_post>Topto.MeanNspk_optoOFF_post);
Opto_post_sess = Topto.OPTO_POST;
Opto_post_Nstim = Topto.OPTO_POST_Nstim;

Ncell_opto =  sum(Topto.H_opto_post)
Ncell_opto_inib = sum(Opto_inib)
Ncell_opto_exct = sum(Opto_exct)

%% From T.ephys
%% TTEST
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

%% From Tephys_z
%% Z-SCORE for Fig2: ALL cor trial 
% rename var from Tephys_z
z_exct_pufAll =Tephys_z.zpuf_exct; 
z_exct_delAll =Tephys_z.zdel_exct;  
z_exct_resAll =Tephys_z.zres_exct;
z_inib_pufAll =Tephys_z.zpuf_inib; 
z_inib_delAll =Tephys_z.zdel_inib;  
z_inib_resAll =Tephys_z.zres_inib;

z_exctAll      =  z_exct_pufAll | z_exct_delAll | z_exct_resAll;   % All cell that are sig Excited in task 
z_inibAll      =  z_inib_pufAll | z_inib_delAll | z_inib_resAll; 
z_exct         =  z_exctAll & ~z_inibAll;  % Excite ONLY 
z_inib         = ~z_exctAll &  z_inibAll; % Inib ONLY 

z_complex      =  z_exctAll &  z_inibAll; % Complex ONLY 
z_nosig        = ~z_exctAll & ~z_inibAll; % NoSig ONLY 

z_exct_UniMod_1puf      =  z_exct &  z_exct_pufAll & ~z_exct_delAll  & ~z_exct_resAll; 
z_exct_UniMod_2del      =  z_exct & ~z_exct_pufAll &  z_exct_delAll  & ~z_exct_resAll; 
z_exct_UniMod_3res      =  z_exct & ~z_exct_pufAll & ~z_exct_delAll  &  z_exct_resAll; 
z_inib_UniMod_1puf      =  z_inib &  z_inib_pufAll & ~z_inib_delAll  & ~z_inib_resAll; 
z_inib_UniMod_2del      =  z_inib & ~z_inib_pufAll &  z_inib_delAll  & ~z_inib_resAll; 
z_inib_UniMod_3res      =  z_inib & ~z_inib_pufAll & ~z_inib_delAll  &  z_inib_resAll; 

z_exct_biMod_12pd      =  z_exct &  z_exct_pufAll &  z_exct_delAll  & ~z_exct_resAll; 
z_exct_biMod_13pr      =  z_exct &  z_exct_pufAll & ~z_exct_delAll  &  z_exct_resAll; 
z_exct_biMod_23dr      =  z_exct & ~z_exct_pufAll &  z_exct_delAll  &  z_exct_resAll; 
z_inib_biMod_12pd      =  z_inib &  z_inib_pufAll &  z_inib_delAll  & ~z_inib_resAll; 
z_inib_biMod_13pr      =  z_inib &  z_inib_pufAll & ~z_inib_delAll  &  z_inib_resAll; 
z_inib_biMod_23dr      =  z_inib & ~z_inib_pufAll &  z_inib_delAll  &  z_inib_resAll; 

z_exct_triMod_123      =  z_exct &  z_exct_pufAll &  z_exct_delAll  &  z_exct_resAll; 
z_inib_triMod_123      =  z_inib &  z_inib_pufAll &  z_inib_delAll  &  z_inib_resAll;

Tcombo_z= [];
Tcombo_z= addvars(Tcoord(:,1:6), VMVL, AreaID,...
Opto_inib, Opto_exct, Opto_post_sess, Opto_post_Nstim, ... 
z_exct_pufAll, z_exct_delAll, z_exct_resAll, ... 
z_inib_pufAll, z_inib_delAll, z_inib_resAll, ...
z_exctAll, z_inibAll, z_exct, z_inib, z_complex, z_nosig, ...
z_exct_UniMod_1puf, z_exct_UniMod_2del, z_exct_UniMod_3res, ... 
z_inib_UniMod_1puf, z_inib_UniMod_2del, z_inib_UniMod_3res, ...
z_exct_biMod_12pd, z_exct_biMod_13pr, z_exct_biMod_23dr, ... 
z_inib_biMod_12pd, z_inib_biMod_13pr, z_inib_biMod_23dr, ...
z_exct_triMod_123, z_inib_triMod_123); 

Tcombo_z(1,:);
save('D:\JC_Analysis\listcell.mat','Tcombo_z', '-append');
disp('Tcombo_z saved')


%% From TcLvR_z
%CONTRA CELLS of each types: exct vs Inib vs Complex vs NotSig 
%% Z-SCORE for Fig2: CONTRA trials 
% rename var from TcLvR_z
zcont_exct_pufAll =TcLvR_z.cont_zpuf_exct; 
zcont_exct_delAll =TcLvR_z.cont_zdel_exct;  
zcont_exct_resAll =TcLvR_z.cont_zres_exct;
zcont_inib_pufAll =TcLvR_z.cont_zpuf_inib; 
zcont_inib_delAll =TcLvR_z.cont_zdel_inib;  
zcont_inib_resAll =TcLvR_z.cont_zres_inib;

zcont_exctAll      =  zcont_exct_pufAll | zcont_exct_delAll | zcont_exct_resAll;   % All cell that are sig Excited in task 
zcont_inibAll      =  zcont_inib_pufAll | zcont_inib_delAll | zcont_inib_resAll; 
zcont_exct         =  zcont_exctAll & ~zcont_inibAll;  % Excite ONLY 
zcont_inib         = ~zcont_exctAll &  zcont_inibAll; % Inib ONLY 

zcont_complex      =  zcont_exctAll &  zcont_inibAll; % Complex ONLY 
zcont_nosig        = ~zcont_exctAll & ~zcont_inibAll; % NoSig ONLY 

zcont_exct_UniMod_1puf      =  zcont_exct &  zcont_exct_pufAll & ~zcont_exct_delAll  & ~zcont_exct_resAll; 
zcont_exct_UniMod_2del      =  zcont_exct & ~zcont_exct_pufAll &  zcont_exct_delAll  & ~zcont_exct_resAll; 
zcont_exct_UniMod_3res      =  zcont_exct & ~zcont_exct_pufAll & ~zcont_exct_delAll  &  zcont_exct_resAll; 
zcont_inib_UniMod_1puf      =  zcont_inib &  zcont_inib_pufAll & ~zcont_inib_delAll  & ~zcont_inib_resAll; 
zcont_inib_UniMod_2del      =  zcont_inib & ~zcont_inib_pufAll &  zcont_inib_delAll  & ~zcont_inib_resAll; 
zcont_inib_UniMod_3res      =  zcont_inib & ~zcont_inib_pufAll & ~zcont_inib_delAll  &  zcont_inib_resAll; 

zcont_exct_biMod_12pd      =  zcont_exct &  zcont_exct_pufAll &  zcont_exct_delAll  & ~zcont_exct_resAll; 
zcont_exct_biMod_13pr      =  zcont_exct &  zcont_exct_pufAll & ~zcont_exct_delAll  &  zcont_exct_resAll; 
zcont_exct_biMod_23dr      =  zcont_exct & ~zcont_exct_pufAll &  zcont_exct_delAll  &  zcont_exct_resAll; 
zcont_inib_biMod_12pd      =  zcont_inib &  zcont_inib_pufAll &  zcont_inib_delAll  & ~zcont_inib_resAll; 
zcont_inib_biMod_13pr      =  zcont_inib &  zcont_inib_pufAll & ~zcont_inib_delAll  &  zcont_inib_resAll; 
zcont_inib_biMod_23dr      =  zcont_inib & ~zcont_inib_pufAll &  zcont_inib_delAll  &  zcont_inib_resAll; 

zcont_exct_triMod_123      =  zcont_exct &  zcont_exct_pufAll &  zcont_exct_delAll  &  zcont_exct_resAll; 
zcont_inib_triMod_123      =  zcont_inib &  zcont_inib_pufAll &  zcont_inib_delAll  &  zcont_inib_resAll;

Tcombo_zcont= [];
Tcombo_zcont= addvars(Tcoord(:,1:6), VMVL, AreaID,...
Opto_inib, Opto_exct, Opto_post_sess, Opto_post_Nstim, ...  
zcont_exct_pufAll, zcont_exct_delAll, zcont_exct_resAll, ... 
zcont_inib_pufAll, zcont_inib_delAll, zcont_inib_resAll, ...
zcont_exctAll, zcont_inibAll, zcont_exct, zcont_inib, zcont_complex, zcont_nosig, ...
zcont_exct_UniMod_1puf, zcont_exct_UniMod_2del, zcont_exct_UniMod_3res, ... 
zcont_inib_UniMod_1puf, zcont_inib_UniMod_2del, zcont_inib_UniMod_3res, ...
zcont_exct_biMod_12pd, zcont_exct_biMod_13pr, zcont_exct_biMod_23dr, ... 
zcont_inib_biMod_12pd, zcont_inib_biMod_13pr, zcont_inib_biMod_23dr, ...
zcont_exct_triMod_123, zcont_inib_triMod_123); 

Tcombo_zcont(1,:);
save('D:\JC_Analysis\listcell.mat','Tcombo_zcont', '-append');
disp('Tcombo_zcont saved')

%% From TcLvR_z
%IPSI CELLS of each types: exct vs Inib vs Complex vs NotSig 
%% Z-SCORE for Fig2: IPSI trials 
% rename var from TcLvR_z
zipsi_exct_pufAll =TcLvR_z.ipsi_zpuf_exct; 
zipsi_exct_delAll =TcLvR_z.ipsi_zdel_exct;  
zipsi_exct_resAll =TcLvR_z.ipsi_zres_exct;
zipsi_inib_pufAll =TcLvR_z.ipsi_zpuf_inib; 
zipsi_inib_delAll =TcLvR_z.ipsi_zdel_inib;  
zipsi_inib_resAll =TcLvR_z.ipsi_zres_inib;

zipsi_exctAll      =  zipsi_exct_pufAll | zipsi_exct_delAll | zipsi_exct_resAll;   % All cell that are sig Excited in task 
zipsi_inibAll      =  zipsi_inib_pufAll | zipsi_inib_delAll | zipsi_inib_resAll; 
zipsi_exct         =  zipsi_exctAll & ~zipsi_inibAll;  % Excite ONLY 
zipsi_inib         = ~zipsi_exctAll &  zipsi_inibAll; % Inib ONLY 

zipsi_complex      =  zipsi_exctAll &  zipsi_inibAll; % Complex ONLY 
zipsi_nosig        = ~zipsi_exctAll & ~zipsi_inibAll; % NoSig ONLY 

zipsi_exct_UniMod_1puf      =  zipsi_exct &  zipsi_exct_pufAll & ~zipsi_exct_delAll  & ~zipsi_exct_resAll; 
zipsi_exct_UniMod_2del      =  zipsi_exct & ~zipsi_exct_pufAll &  zipsi_exct_delAll  & ~zipsi_exct_resAll; 
zipsi_exct_UniMod_3res      =  zipsi_exct & ~zipsi_exct_pufAll & ~zipsi_exct_delAll  &  zipsi_exct_resAll; 
zipsi_inib_UniMod_1puf      =  zipsi_inib &  zipsi_inib_pufAll & ~zipsi_inib_delAll  & ~zipsi_inib_resAll; 
zipsi_inib_UniMod_2del      =  zipsi_inib & ~zipsi_inib_pufAll &  zipsi_inib_delAll  & ~zipsi_inib_resAll; 
zipsi_inib_UniMod_3res      =  zipsi_inib & ~zipsi_inib_pufAll & ~zipsi_inib_delAll  &  zipsi_inib_resAll; 

zipsi_exct_biMod_12pd      =  zipsi_exct &  zipsi_exct_pufAll &  zipsi_exct_delAll  & ~zipsi_exct_resAll; 
zipsi_exct_biMod_13pr      =  zipsi_exct &  zipsi_exct_pufAll & ~zipsi_exct_delAll  &  zipsi_exct_resAll; 
zipsi_exct_biMod_23dr      =  zipsi_exct & ~zipsi_exct_pufAll &  zipsi_exct_delAll  &  zipsi_exct_resAll; 
zipsi_inib_biMod_12pd      =  zipsi_inib &  zipsi_inib_pufAll &  zipsi_inib_delAll  & ~zipsi_inib_resAll; 
zipsi_inib_biMod_13pr      =  zipsi_inib &  zipsi_inib_pufAll & ~zipsi_inib_delAll  &  zipsi_inib_resAll; 
zipsi_inib_biMod_23dr      =  zipsi_inib & ~zipsi_inib_pufAll &  zipsi_inib_delAll  &  zipsi_inib_resAll; 

zipsi_exct_triMod_123      =  zipsi_exct &  zipsi_exct_pufAll &  zipsi_exct_delAll  &  zipsi_exct_resAll; 
zipsi_inib_triMod_123      =  zipsi_inib &  zipsi_inib_pufAll &  zipsi_inib_delAll  &  zipsi_inib_resAll;

Tcombo_zipsi= [];
Tcombo_zipsi= addvars(Tcoord(:,1:6), VMVL, AreaID, ...
Opto_inib, Opto_exct, Opto_post_sess, Opto_post_Nstim, ...   
zipsi_exct_pufAll, zipsi_exct_delAll, zipsi_exct_resAll, ... 
zipsi_inib_pufAll, zipsi_inib_delAll, zipsi_inib_resAll, ...
zipsi_exctAll, zipsi_inibAll, zipsi_exct, zipsi_inib, zipsi_complex, zipsi_nosig, ...
zipsi_exct_UniMod_1puf, zipsi_exct_UniMod_2del, zipsi_exct_UniMod_3res, ... 
zipsi_inib_UniMod_1puf, zipsi_inib_UniMod_2del, zipsi_inib_UniMod_3res, ...
zipsi_exct_biMod_12pd, zipsi_exct_biMod_13pr, zipsi_exct_biMod_23dr, ... 
zipsi_inib_biMod_12pd, zipsi_inib_biMod_13pr, zipsi_inib_biMod_23dr, ...
zipsi_exct_triMod_123, zipsi_inib_triMod_123); 

Tcombo_zipsi(1,:);
save('D:\JC_Analysis\listcell.mat','Tcombo_zipsi', '-append');
disp('Tcombo_zipsi saved')


