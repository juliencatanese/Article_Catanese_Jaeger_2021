function [sort_Mat2sort] = SortPeakSDF_JCfun(Mat2sort, sort_type, parfig)
% function [sort_SNMA, sort_SNSA, sort_SNSemA] = SortPeakSDF_JCfun(SNMA, SNSA, SNSemA, sort_type)
% sort_type = 'ascend' or 'descend'
XsortEpoch=parfig.XsortEpoch; 
clear IDX IDXnan sort_Mat2sort;
IDXnan = zeros(size(Mat2sort,1),1);
for irow= 1:size(Mat2sort,1);
    if ~isnan(Mat2sort(irow,:));
        IDX(irow) =  min(find(abs(Mat2sort(irow, XsortEpoch))==max(abs(Mat2sort(irow, [XsortEpoch])))));
    else
        IDXnan(irow)=1;
    end
end

Mat2sort(find(IDXnan),:) =[];
IDX(find(IDXnan)) =[]; 

Mat2sort = [IDX' Mat2sort];

sort_Mat2sort = sortrows(Mat2sort,sort_type);

sort_Mat2sort = sort_Mat2sort(2:end,:);
