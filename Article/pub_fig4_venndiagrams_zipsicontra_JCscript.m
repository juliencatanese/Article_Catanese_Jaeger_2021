

clear all, close all, 

load('listcell.mat')

Nzcont_exct_pufAll = sum(Tcombo.VMVL  & Tcombo_zcont.zcont_exct_pufAll); 
Nzcont_exct_delAll = sum(Tcombo.VMVL  & Tcombo_zcont.zcont_exct_delAll);  
Nzcont_exct_resAll = sum(Tcombo.VMVL  & Tcombo_zcont.zcont_exct_resAll);
Nzcont_inib_pufAll = sum(Tcombo.VMVL  & Tcombo_zcont.zcont_inib_pufAll); 
Nzcont_inib_delAll = sum(Tcombo.VMVL  & Tcombo_zcont.zcont_inib_delAll);  
Nzcont_inib_resAll = sum(Tcombo.VMVL  & Tcombo_zcont.zcont_inib_resAll);

Nzcont_exctAll      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_exctAll)   % All cell that are sig Excited in task 
Nzcont_inibAll      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_inibAll) 
Nzcont_exct         =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_exct)  % Excite ONLY 
Nzcont_inib         =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_inib) % Inib ONLY 
Nzcont_complex      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_complex) % Complex ONLY 
Nzcont_nosig        =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_nosig)  % NoSig ONLY 

Nzcont_exct_UniMod_1puf      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_exct_UniMod_1puf)
Nzcont_exct_UniMod_2del      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_exct_UniMod_2del)
Nzcont_exct_UniMod_3res      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_exct_UniMod_3res)
Nzcont_inib_UniMod_1puf      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_inib_UniMod_1puf )
Nzcont_inib_UniMod_2del      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_inib_UniMod_2del)
Nzcont_inib_UniMod_3res      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_inib_UniMod_3res )

Nzcont_exct_biMod_12pd      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_exct_biMod_12pd)
Nzcont_exct_biMod_13pr      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_exct_biMod_13pr)
Nzcont_exct_biMod_23dr      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_exct_biMod_23dr)
Nzcont_inib_biMod_12pd      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_inib_biMod_12pd)
Nzcont_inib_biMod_13pr      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_inib_biMod_13pr)
Nzcont_inib_biMod_23dr      =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_inib_biMod_23dr)

Nzcont_exct_triMod_123    =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_exct_triMod_123)
Nzcont_inib_triMod_123    =  sum(Tcombo.VMVL  & Tcombo_zcont.zcont_inib_triMod_123)



Nzipsi_exct_pufAll = sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_exct_pufAll); 
Nzipsi_exct_delAll = sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_exct_delAll);  
Nzipsi_exct_resAll = sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_exct_resAll);
Nzipsi_inib_pufAll = sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_inib_pufAll); 
Nzipsi_inib_delAll = sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_inib_delAll);  
Nzipsi_inib_resAll = sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_inib_resAll);

Nzipsi_exctAll      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_exctAll)   % All cell that are sig Excited in task 
Nzipsi_inibAll      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_inibAll) 
Nzipsi_exct         =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_exct)  % Excite ONLY 
Nzipsi_inib         =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_inib) % Inib ONLY 
Nzipsi_complex      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_complex) % Complex ONLY 
Nzipsi_nosig        =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_nosig)  % NoSig ONLY 

Nzipsi_exct_UniMod_1puf      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_exct_UniMod_1puf)
Nzipsi_exct_UniMod_2del      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_exct_UniMod_2del)
Nzipsi_exct_UniMod_3res      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_exct_UniMod_3res)
Nzipsi_inib_UniMod_1puf      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_inib_UniMod_1puf )
Nzipsi_inib_UniMod_2del      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_inib_UniMod_2del)
Nzipsi_inib_UniMod_3res      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_inib_UniMod_3res )

Nzipsi_exct_biMod_12pd      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_exct_biMod_12pd)
Nzipsi_exct_biMod_13pr      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_exct_biMod_13pr)
Nzipsi_exct_biMod_23dr      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_exct_biMod_23dr)
Nzipsi_inib_biMod_12pd      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_inib_biMod_12pd)
Nzipsi_inib_biMod_13pr      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_inib_biMod_13pr)
Nzipsi_inib_biMod_23dr      =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_inib_biMod_23dr)

Nzipsi_exct_triMod_123    =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_exct_triMod_123)
Nzipsi_inib_triMod_123    =  sum(Tcombo.VMVL  & Tcombo_zipsi.zipsi_inib_triMod_123)

%% PLOT vennDiagrams (VennX)  
% data = [1 12 2 23 3 13 123]

% IPSI exct v inib 
data= [Nzipsi_exct, Nzipsi_complex, Nzipsi_inib]
resolution= 0.05
colmap='bone'
vennX(data, resolution, colmap)
title('ipsi')
colormapeditor

% CONTRA exct v inib 
data= [Nzcont_exct, Nzcont_complex, Nzcont_inib]
resolution= 0.05
colmap='bone'
vennX(data, resolution, colmap)
title('contra')
colormapeditor

% IPSI Types
data= [Nzipsi_exct_UniMod_1puf, Nzipsi_exct_biMod_12pd, Nzipsi_exct_UniMod_2del, Nzipsi_exct_biMod_23dr, Nzipsi_exct_UniMod_3res, Nzipsi_exct_biMod_13pr, Nzipsi_exct_triMod_123 ] 
resolution= 0.05
colmap='bone'
vennX(data, resolution, colmap)
title('ipsi')
colormapeditor

% CONTRA Types
data= [Nzcont_exct_UniMod_1puf, Nzcont_exct_biMod_12pd, Nzcont_exct_UniMod_2del, Nzcont_exct_biMod_23dr, Nzcont_exct_UniMod_3res, Nzcont_exct_biMod_13pr, Nzcont_exct_triMod_123 ] 
resolution= 0.05
colmap='bone'
vennX(data, resolution, colmap)
title('contra')
colormapeditor