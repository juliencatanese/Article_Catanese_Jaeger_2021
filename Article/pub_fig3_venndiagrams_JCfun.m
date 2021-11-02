function pub_fig3_venndiagrams_JCfun(Tzside, side, colmap, resolution)
% function pub_fig3_venndiagrams_JCfun(Tzside, colmap, resolution)
% plot venndiagrams
% colmap = colormap set (e.g. 'jet', 'bone', 'autumn')
% written by julien Catanese 3/31/2019

if side == 'cont'
    Nz_exct_pufAll = sum(Tzside.VMVL  & Tzside.zcont_exct_pufAll);
    Nz_exct_delAll = sum(Tzside.VMVL  & Tzside.zcont_exct_delAll);
    Nz_exct_resAll = sum(Tzside.VMVL  & Tzside.zcont_exct_resAll);
    Nz_inib_pufAll = sum(Tzside.VMVL  & Tzside.zcont_inib_pufAll);
    Nz_inib_delAll = sum(Tzside.VMVL  & Tzside.zcont_inib_delAll);
    Nz_inib_resAll = sum(Tzside.VMVL  & Tzside.zcont_inib_resAll);
    
    Nz_exctAll      =  sum(Tzside.VMVL  & Tzside.zcont_exctAll)   % All cell that are sig Excited in task
    Nz_inibAll      =  sum(Tzside.VMVL  & Tzside.zcont_inibAll)
    Nz_exct         =  sum(Tzside.VMVL  & Tzside.zcont_exct)  % Excite ONLY
    Nz_inib         =  sum(Tzside.VMVL  & Tzside.zcont_inib) % Inib ONLY
    Nz_complex      =  sum(Tzside.VMVL  & Tzside.zcont_complex) % Complex ONLY
    Nz_nosig        =  sum(Tzside.VMVL  & Tzside.zcont_nosig)  % NoSig ONLY
    
    Nz_exct_UniMod_1puf      =  sum(Tzside.VMVL  & Tzside.zcont_exct_UniMod_1puf)
    Nz_exct_UniMod_2del      =  sum(Tzside.VMVL  & Tzside.zcont_exct_UniMod_2del)
    Nz_exct_UniMod_3res      =  sum(Tzside.VMVL  & Tzside.zcont_exct_UniMod_3res)
    Nz_inib_UniMod_1puf      =  sum(Tzside.VMVL  & Tzside.zcont_inib_UniMod_1puf )
    Nz_inib_UniMod_2del      =  sum(Tzside.VMVL  & Tzside.zcont_inib_UniMod_2del)
    Nz_inib_UniMod_3res      =  sum(Tzside.VMVL  & Tzside.zcont_inib_UniMod_3res )
    
    Nz_exct_biMod_12pd      =  sum(Tzside.VMVL  & Tzside.zcont_exct_biMod_12pd)
    Nz_exct_biMod_13pr      =  sum(Tzside.VMVL  & Tzside.zcont_exct_biMod_13pr)
    Nz_exct_biMod_23dr      =  sum(Tzside.VMVL  & Tzside.zcont_exct_biMod_23dr)
    Nz_inib_biMod_12pd      =  sum(Tzside.VMVL  & Tzside.zcont_inib_biMod_12pd)
    Nz_inib_biMod_13pr      =  sum(Tzside.VMVL  & Tzside.zcont_inib_biMod_13pr)
    Nz_inib_biMod_23dr      =  sum(Tzside.VMVL  & Tzside.zcont_inib_biMod_23dr)
    
    Nz_exct_triMod_123    =  sum(Tzside.VMVL  & Tzside.zcont_exct_triMod_123)
    Nz_inib_triMod_123    =  sum(Tzside.VMVL  & Tzside.zcont_inib_triMod_123)
    
    
elseif side == 'ipsi'
    Nz_exct_pufAll = sum(Tzside.VMVL  & Tzside.zipsi_exct_pufAll);
    Nz_exct_delAll = sum(Tzside.VMVL  & Tzside.zipsi_exct_delAll);
    Nz_exct_resAll = sum(Tzside.VMVL  & Tzside.zipsi_exct_resAll);
    Nz_inib_pufAll = sum(Tzside.VMVL  & Tzside.zipsi_inib_pufAll);
    Nz_inib_delAll = sum(Tzside.VMVL  & Tzside.zipsi_inib_delAll);
    Nz_inib_resAll = sum(Tzside.VMVL  & Tzside.zipsi_inib_resAll);
    
    Nz_exctAll      =  sum(Tzside.VMVL  & Tzside.zipsi_exctAll)   % All cell that are sig Excited in task
    Nz_inibAll      =  sum(Tzside.VMVL  & Tzside.zipsi_inibAll)
    Nz_exct         =  sum(Tzside.VMVL  & Tzside.zipsi_exct)  % Excite ONLY
    Nz_inib         =  sum(Tzside.VMVL  & Tzside.zipsi_inib) % Inib ONLY
    Nz_complex      =  sum(Tzside.VMVL  & Tzside.zipsi_complex) % Complex ONLY
    Nz_nosig        =  sum(Tzside.VMVL  & Tzside.zipsi_nosig)  % NoSig ONLY
    
    Nz_exct_UniMod_1puf      =  sum(Tzside.VMVL  & Tzside.zipsi_exct_UniMod_1puf)
    Nz_exct_UniMod_2del      =  sum(Tzside.VMVL  & Tzside.zipsi_exct_UniMod_2del)
    Nz_exct_UniMod_3res      =  sum(Tzside.VMVL  & Tzside.zipsi_exct_UniMod_3res)
    Nz_inib_UniMod_1puf      =  sum(Tzside.VMVL  & Tzside.zipsi_inib_UniMod_1puf )
    Nz_inib_UniMod_2del      =  sum(Tzside.VMVL  & Tzside.zipsi_inib_UniMod_2del)
    Nz_inib_UniMod_3res      =  sum(Tzside.VMVL  & Tzside.zipsi_inib_UniMod_3res )
    
    Nz_exct_biMod_12pd      =  sum(Tzside.VMVL  & Tzside.zipsi_exct_biMod_12pd)
    Nz_exct_biMod_13pr      =  sum(Tzside.VMVL  & Tzside.zipsi_exct_biMod_13pr)
    Nz_exct_biMod_23dr      =  sum(Tzside.VMVL  & Tzside.zipsi_exct_biMod_23dr)
    Nz_inib_biMod_12pd      =  sum(Tzside.VMVL  & Tzside.zipsi_inib_biMod_12pd)
    Nz_inib_biMod_13pr      =  sum(Tzside.VMVL  & Tzside.zipsi_inib_biMod_13pr)
    Nz_inib_biMod_23dr      =  sum(Tzside.VMVL  & Tzside.zipsi_inib_biMod_23dr)
    
    Nz_exct_triMod_123    =  sum(Tzside.VMVL  & Tzside.zipsi_exct_triMod_123)
    Nz_inib_triMod_123    =  sum(Tzside.VMVL  & Tzside.zipsi_inib_triMod_123)
end

%% plot VENNDIAGRAM
% VennX: data = [1 12 2 23 3 13 123]

% +/-/complex/
data= [Nz_exct, Nz_complex, Nz_inib]
vennX(data, resolution, colmap)
colormapeditor

data= [Nz_exct_UniMod_1puf, Nz_exct_biMod_12pd, Nz_exct_UniMod_2del, Nz_exct_biMod_23dr, Nz_exct_UniMod_3res, Nz_exct_biMod_13pr, Nz_exct_triMod_123 ]
vennX(data, resolution, colmap)
colormapeditor