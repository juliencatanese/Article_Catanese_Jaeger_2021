
function [Sig_idx, BarSig] = Select_sigCell_1std_JCfun_v2(sdf_mean, sdf_std, Sig_idx, BarSig, pre, post, ncell)

Baselin_mean_Av = mean(sdf_mean(163:2250));
Baselin_Std_Av = mean(sdf_std(163:2250));
thr = (Baselin_mean_Av + Baselin_Std_Av);
Cell = sdf_mean;

idxB = find(Cell >= thr);
BarBin= zeros(1,pre+post+1);
BarBin(idxB)= 1;

if  ~isempty(idxB)  
    BarSig = [ BarSig ; BarBin];
    Sig_idx = [Sig_idx ncell] ; % ncell is the idx in SMA.  
end