function pub_fig_venndiagrams_zcor_JCfun(Tcombo, colmap)
% function pub_fig_venndiagrams_zcor_JCfun(Tcombo, colmap)
% plot venndiagrams 
% colmap = colormap set (e.g. 'jet', 'bone', 'autumn')
% written by julien Catanese 3/31/2019 

Nz_exct_pufAll = sum(Tcombo.VMVL  & Tcombo_z.z_exct_pufAll); 
Nz_exct_delAll = sum(Tcombo.VMVL  & Tcombo_z.z_exct_delAll);  
Nz_exct_resAll = sum(Tcombo.VMVL  & Tcombo_z.z_exct_resAll);
Nz_inib_pufAll = sum(Tcombo.VMVL  & Tcombo_z.z_inib_pufAll); 
Nz_inib_delAll = sum(Tcombo.VMVL  & Tcombo_z.z_inib_delAll);  
Nz_inib_resAll = sum(Tcombo.VMVL  & Tcombo_z.z_inib_resAll);

Nz_exctAll      =  sum(Tcombo.VMVL  & Tcombo_z.z_exctAll)   % All cell that are sig Excited in task 
Nz_inibAll      =  sum(Tcombo.VMVL  & Tcombo_z.z_inibAll) 
Nz_exct         =  sum(Tcombo.VMVL  & Tcombo_z.z_exct)  % Excite ONLY 
Nz_inib         =  sum(Tcombo.VMVL  & Tcombo_z.z_inib) % Inib ONLY 
Nz_complex      =  sum(Tcombo.VMVL  & Tcombo_z.z_complex) % Complex ONLY 
Nz_nosig        =  sum(Tcombo.VMVL  & Tcombo_z.z_nosig)  % NoSig ONLY 

Nz_exct_UniMod_1puf      =  sum(Tcombo.VMVL  & Tcombo_z.z_exct_UniMod_1puf)
Nz_exct_UniMod_2del      =  sum(Tcombo.VMVL  & Tcombo_z.z_exct_UniMod_2del)
Nz_exct_UniMod_3res      =  sum(Tcombo.VMVL  & Tcombo_z.z_exct_UniMod_3res)
Nz_inib_UniMod_1puf      =  sum(Tcombo.VMVL  & Tcombo_z.z_inib_UniMod_1puf )
Nz_inib_UniMod_2del      =  sum(Tcombo.VMVL  & Tcombo_z.z_inib_UniMod_2del)
Nz_inib_UniMod_3res      =  sum(Tcombo.VMVL  & Tcombo_z.z_inib_UniMod_3res )

Nz_exct_biMod_12pd      =  sum(Tcombo.VMVL  & Tcombo_z.z_exct_biMod_12pd)
Nz_exct_biMod_13pr      =  sum(Tcombo.VMVL  & Tcombo_z.z_exct_biMod_13pr)
Nz_exct_biMod_23dr      =  sum(Tcombo.VMVL  & Tcombo_z.z_exct_biMod_23dr)
Nz_inib_biMod_12pd      =  sum(Tcombo.VMVL  & Tcombo_z.z_inib_biMod_12pd)
Nz_inib_biMod_13pr      =  sum(Tcombo.VMVL  & Tcombo_z.z_inib_biMod_13pr)
Nz_inib_biMod_23dr      =  sum(Tcombo.VMVL  & Tcombo_z.z_inib_biMod_23dr)

Nz_exct_triMod_123    =  sum(Tcombo.VMVL  & Tcombo_z.z_exct_triMod_123)
Nz_inib_triMod_123    =  sum(Tcombo.VMVL  & Tcombo_z.z_inib_triMod_123)

%% plot VENNDIAGRAM
% VennX: data = [1 12 2 23 3 13 123]

% +/-/complex/
data= [Nz_exct, Nz_complex, Nz_inib]
resolution= 0.05
colmap='bone'
vennX(data, resolution, colmap)
colormapeditor

data= [Nz_exct_UniMod_1puf, Nz_exct_biMod_12pd, Nz_exct_UniMod_2del, Nz_exct_biMod_23dr, Nz_exct_UniMod_3res, Nz_exct_biMod_13pr, Nz_exct_triMod_123 ] 
resolution= 0.05
colmap='bone'
vennX(data, resolution, colmap)
colormapeditor