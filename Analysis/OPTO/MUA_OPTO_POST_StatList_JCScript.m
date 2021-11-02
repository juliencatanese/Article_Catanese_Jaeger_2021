% MUA_OPTO_POST_StatList_JCScript
% Loop through dataset and Call function: 
% [HOFFon,POFFon, HONoff, PONoff]= MUA_OPTO_POST_ttest_JCfun(ChanID, dispFig)
% Compare Nspike>thr in OPTO vs NON-OPTO timeStim during the POSTtask epoch 
% Generate a list of ranktest results for each Channel of each Session
% [SessID_all ChanID_all HOFFon_all POFFon_all]
% Written by JC 11/13/2018
% Last modified by JC 2/5/2019

clear all

cd('D:\JC_Analysis');
SessList = dir(['*JCVGAT1*/*post*']);
NSess= max(size(SessList)) % Number of Sessions
recap=1
SessID_all=[]; ChanID_all=[]; HOFFon_all=[]; POFFon_all=[];  HONoff_all=[];  PONoff_all=[];


%% loop trhough all Sessions named "taskopto"
for ns=recap:NSess;    
    SessID_full= SessList(ns).name;
    SessPath=SessList(ns).folder;
    cd([SessPath '\' SessID_full]);
    SessID=SessID_full(1:12)
    ChID_list = dir(['S*Ch*_sub.mat']);
    for nc= 1:max(size(ChID_list))
        ChanID=ChID_list(nc).name(1:5); 
        dispFig=0; 
        [HOFFon,POFFon, HONoff, PONoff]= MUA_OPTO_POST_ttest_JCfun(ChanID, dispFig); disp(['HOFFon = ' num2str(HOFFon) ' for '  ChanID ]); 
        SessID_all=[SessID_all; SessID];
        ChanID_all=[ChanID_all; ChanID];
        HOFFon_all=[HOFFon_all; HOFFon];
        POFFon_all=[POFFon_all; POFFon];
        HONoff_all=[HONoff_all; HONoff];
        PONoff_all=[PONoff_all; PONoff];
        
    end
    
end

%%
Session = cellstr(SessID_all); 
AA= table(Session , ChanID_all, HOFFon_all, HONoff_all, POFFon_all, PONoff_all )

ListSess=unique(AA.Session)
for ii=1:size(ListSess,1)
AA(strcmp(AA.Session, ListSess{ii}),:)
[sum(AA.HOFFon_all(strcmp(AA.Session, ListSess{ii}),:)) sum(AA.HONoff_all(strcmp(AA.Session, ListSess{ii}),:))]
end
