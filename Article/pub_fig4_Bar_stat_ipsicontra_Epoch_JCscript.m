
% pub_fig4_Bar_stat_ipsicontra_Epoch_JCscript 
% contra_dz = index of contra cell including comlpex 
% last modified 11/17/2019 by Catanese Julien 

close all, figure,


data_mean =[abs(mean(Tfig3_RvL.dzSMA_MAX_puf(contra_dz))), abs(mean(Tfig3_RvL.dzSMA_MIN_puf(ipsi_dz))); ...
    abs(mean(Tfig3_RvL.dzSMA_MAX_del(contra_dz))), abs(mean(Tfig3_RvL.dzSMA_MIN_del(ipsi_dz))); ...
    abs(mean(Tfig3_RvL.dzSMA_MAX_res(contra_dz))), abs(mean(Tfig3_RvL.dzSMA_MIN_res(ipsi_dz)))];
data_sem = [abs(std(Tfig3_RvL.dzSMA_MAX_puf(contra_dz)/sqrt(size(contra_dz,1)))), abs(std(Tfig3_RvL.dzSMA_MIN_puf(ipsi_dz)/sqrt(size(ipsi_dz,1)))); ...
    abs(std(Tfig3_RvL.dzSMA_MAX_del(contra_dz)/sqrt(size(contra_dz,1)))), abs(std(Tfig3_RvL.dzSMA_MIN_del(ipsi_dz)/sqrt(size(ipsi_dz,1))));...
    abs(std(Tfig3_RvL.dzSMA_MAX_res(contra_dz)/sqrt(size(contra_dz,1)))), abs(std(Tfig3_RvL.dzSMA_MIN_res(ipsi_dz)/sqrt(size(ipsi_dz,1))))];

% non-param test
[Ppuf Hpuf] = ranksum(Tfig3_RvL.dzSMA_MAX_puf(contra_dz), -Tfig3_RvL.dzSMA_MIN_puf(ipsi_dz))
[Pdel Hdel] = ranksum(Tfig3_RvL.dzSMA_MAX_del(contra_dz), -Tfig3_RvL.dzSMA_MIN_del(ipsi_dz))
[Pres Hres] = ranksum(Tfig3_RvL.dzSMA_MAX_res(contra_dz), -Tfig3_RvL.dzSMA_MIN_res(ipsi_dz))

% param paired ttest
[Hpuf Ppuf] = ttest2(Tfig3_RvL.dzSMA_MAX_puf(contra_dz), -Tfig3_RvL.dzSMA_MIN_puf(ipsi_dz))
[Hdel Pdel] = ttest2(Tfig3_RvL.dzSMA_MAX_del(contra_dz), -Tfig3_RvL.dzSMA_MIN_del(ipsi_dz))
[Hres Pres] = ttest2(Tfig3_RvL.dzSMA_MAX_res(contra_dz), -Tfig3_RvL.dzSMA_MIN_res(ipsi_dz))

bar(data_mean)
hold on,
errorbar([1-0.15 1+0.15; 2-0.15 2+0.15; 3-0.15 3+0.15] ,data_mean, data_sem*3, '.')

