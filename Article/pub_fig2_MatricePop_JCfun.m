function pub_fig2_MatricePop_JCfun(ZorTT, szSMA, Tcombo, parfig)
% function pub_fig2_MatricePop_JCfun(ZorTT, zSMA, parfig)
% OUTPUTS; plot a Matrices (1 colored line per cells) 
% INPUTS: ZorTT choose your way to select cells (Zscore or ttest)
% zSMA or SMA (SdF Mean All) and Tcombo_z or Tcombo (Table of cells)
% parfig: parameters to use for figures. 
% written by Julien Catanese 3/31/2019

if ZorTT=='z_thr'
    parfig.Zthr= 2;
    NszSMA= szSMA./max(abs(szSMA'))';
    nsc=size(szSMA,1);
    
elseif ZorTT=='ttest'
    %  NOT DONE YET
    parfig.title = [SortBY 'Sorted (#Cell='  num2str(nsc) ' Selected by ttest >' num2str(parfig.Zthr) ' trials-' parfig.trial_type{1} ];
    szSMA = zSMA(logical(Tcombo.VMVL) & (type5_Exc | type6_Exc | type7_Exc) ,:);
    szSMA = zSMA(logical(Tcombo.VMVL) & (puf_Exc | del_Exc | res_Exc) ,:);
    NszSMA= szSMA./max(abs(szSMA'))';
    nsc=size(szSMA,1);
end

% plot Matrices/distrib
% parfig.colormap = 'jet';
parfig.sort_direction = 'ascend' ;%parfig.sort_direction = 'descend'
parfig.xlim = [parfig.pre-2250 parfig.pre+750];
parfig.title = [parfig.sort_variable 'Sorted (#Cell='  num2str(nsc) ' Selected by Zthr >' num2str(parfig.Zthr) ' trials-' parfig.trial_type{1} ];

if parfig.sort_variable == 'peak';
    Mat2Sort = NszSMA;
%     parfig.caxis = [0 1];
%     parfig.sort_Xepoch = [parfig.pre-1500 : parfig.pre+750];
elseif parfig.sort_variable == 'Z-tr';
    Mat2Sort = szSMA;
%     parfig.caxis = [-1 12];
%     parfig.sort_Xepoch = [parfig.pre-750 : parfig.pre+1000];
end

[sort_Mat2sort]=pub_comp_SortSDFMat_JCfun(Mat2Sort, parfig);
pub_fig_SeqMat_JCfun(sort_Mat2sort, parfig);





