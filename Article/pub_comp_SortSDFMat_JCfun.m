function [sort_Mat2sort] = pub_comp_SortSDFMat_JCfun(Mat2sort, parfig)
% function [sort_Mat2sort] = pub_comp_SortSDFMat_JCfun(Mat2sort, parfig)
% Necessary to plot sorted Matrice (pub_fig_SeqMat_JCfun.m)
% sort_variable = 'peak' or 'Zthr'
% sort_direction = 'ascend' or 'descend'
% sort_Xepoch =  e.g.  [pre-1500:pre+post-150];
% Zthr= e.g  3 or 10 zscore
% written by Julien Catanese 13MAr2019
%%
sort_variable = parfig.sort_variable;
sort_Xepoch = parfig.sort_Xepoch;
sort_direction = parfig.sort_direction ;
Zthr= parfig.Zthr;

clear IDX IDXnan sort_Mat2sort
IDXnan = zeros(size(Mat2sort,1),1);
for irow= 1:size(Mat2sort,1);
    if ~isnan(Mat2sort(irow,:));
        if sort_variable == 'peak'
            IDX(irow) =  min(find(abs(Mat2sort(irow, sort_Xepoch))==max(abs(Mat2sort(irow, [sort_Xepoch])))));
        elseif sort_variable == 'Z-tr'
            try
                IDX(irow) =  min(find(abs(Mat2sort(irow, sort_Xepoch))>(Zthr+2)))
            catch
                try
                    IDX(irow) =  min(find(abs(Mat2sort(irow, sort_Xepoch))>(Zthr+1)))
                catch
                    try
                        IDX(irow) =  min(find(abs(Mat2sort(irow, sort_Xepoch))>(Zthr)))
                    catch
                        disp('ERROR /2')
                        IDX(irow) =  min(find(abs(Mat2sort(irow, sort_Xepoch))>(Zthr-1)))
                        irow
                    end
                end
            end
        else
            IDXnan(irow)=1;
        end
    end
end
Mat2sort(find(IDXnan),:) =[];
IDX(find(IDXnan)) =[];

Mat2sort = [IDX' Mat2sort];

sort_Mat2sort = sortrows(Mat2sort, sort_direction);

sort_Mat2sort = sort_Mat2sort(2:end,:);

