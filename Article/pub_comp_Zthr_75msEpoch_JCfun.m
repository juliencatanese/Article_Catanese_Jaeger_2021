function  [zbool zvalu] = pub_subT_TaskMod_ZmaxThrEpoch_JCfun(Zthr, zSMA, parfig)
% function pub_subT_TaskMod_ZmaxThrEpoch_JCfun(Zthr, zSMA, parfig) 
% written by Julien Catanese 03/22/2019

GO = parfig.pre; 

% Get the peack (up or Down) for each epoch of the task
zvalu.puf.max = max(zSMA(:,GO-1500:GO-750)')';
zvalu.del.max = max(zSMA(:,GO-750:GO)')';
zvalu.res.max = max(zSMA(:,GO:GO+750)')';
zvalu.puf.min = min(zSMA(:,GO-1500:GO-750)')';
zvalu.del.min = min(zSMA(:,GO-750:GO)')';
zvalu.res.min = min(zSMA(:,GO:GO+750)')';
zvalu.all.abs = max(abs(zSMA(:,200:end-200))')';

Sigexct_lenght_puf =  sum(zSMA(:,GO-1500:GO-750)>  Zthr-1, 2)
Sigexct_lenght_del =  sum(zSMA(:,GO-750:GO)     >  Zthr-1, 2)
Sigexct_lenght_res =  sum(zSMA(:,GO:GO+750)     >  Zthr-1, 2)
Siginib_lenght_puf =  sum(zSMA(:,GO-1500:GO-750)< -Zthr+1, 2)
Siginib_lenght_del =  sum(zSMA(:,GO-750:GO)     < -Zthr+1, 2)
Siginib_lenght_res =  sum(zSMA(:,GO:GO+750)     < -Zthr+1, 2)


zbool.puf.exct = zvalu.puf.max > Zthr  & Sigexct_lenght_puf > 75;
zbool.del.exct = zvalu.del.max > Zthr  & Sigexct_lenght_del > 75;
zbool.res.exct = zvalu.res.max > Zthr  & Sigexct_lenght_res > 75;

zbool.puf.inib = zvalu.puf.min < -Zthr & Siginib_lenght_puf > 75;
zbool.del.inib = zvalu.del.min < -Zthr & Siginib_lenght_del > 75;
zbool.res.inib = zvalu.res.min < -Zthr & Siginib_lenght_res > 75;

% VAR to ADD to TABLE
% load ('listcell.mat'); T1=[]; T2= []; T1=Tephys; Tephys_z = []; 
% T2 = addvars(T1, zSMA_AbsPeak_all, zSMA_MAX_puf, zSMA_MAX_del, zSMA_MAX_res, zSMA_MIN_puf, zSMA_MIN_del, zSMA_MIN_res, zpuf_exct, zdel_exct, zres_exct, zpuf_inib, zdel_inib, zres_inib); 
% Tephys_z = T2;  
% save('D:\JC_Analysis\listcell.mat','Tephys_z', '-append'); 
% Tephys_z(150:155,:) 
% disp('Tephys_z SAVED');


