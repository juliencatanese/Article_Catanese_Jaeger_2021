% pub_table1_get_listcell_JCscript
% create and save a table ('listcell') that contain all cells informations
% ncell, MouseID, Day, DepthZ, ChanID, CluFile2load, SessPath, SessID.  
% save/load('D:\JC_Analysis\listcell.mat','listcell')
% written by Julien Catanese 12/05/2018
% last updated JC 2/7/2019. 

ncell=0; ncell_all = [] ; MouseID_all = []; Day_all =[]; 
SessID_all = []; CluFile2load_all = [] ; SessPath_all = []; 
DepthZ_all = []; ChanID_all = []; CLUST_all = []; clu_str_all = []; 
OptoPost_all = [];

SessList = dir(['*/vgat1*task*']);
NSess= max(size(SessList));

for ns=1:NSess
    SessID= {SessList(ns).name};
    SessPath=SessList(ns).folder;

    cd([SessPath '\' SessID{1}]);
    load('info.mat');
    MouseID = info.info_notes.MouseID;
    Day=info.info_notes.Day;
    DepthZ = info.info_notes.Depth;
    S=SessID{1}; 
    OptoPostID = S(32:39) ;
    
    clust_file = dir('times_*S*Ch*_sub.mat'); ncluf=max(size(clust_file));
    for nchan = 1:ncluf %  channel loop
        
        ChanID=clust_file(nchan).name(7:11);
        load(clust_file(nchan).name);
        Nclust = max(cluster_class(:,1));
        for CLUST= 1:Nclust %Cluster loop within each channel      
            ncell=ncell+1; 
            
            if CLUST>9
                clu_str = ['clu#' num2str(CLUST)]; 
            else 
                clu_str = ['clu#0' num2str(CLUST)];
            end
            CluFile2load =  [clust_file(nchan).name];
            
            
            if OptoPostID == 'nonepost'
                OptoPost = 0 ;
            elseif OptoPostID == 'optopost'
                OptoPost= 1 ;
            end
            
            disp([MouseID ' ' Day ' ' DepthZ ' ' ChanID ' clu#' num2str(CLUST)])
         
            ncell_all = [ncell_all; ncell];            
            MouseID_all = [MouseID_all; MouseID];
            Day_all =[Day_all; Day];
            DepthZ_all = [DepthZ_all; DepthZ];
            ChanID_all = [ChanID_all; ChanID];
            CLUST_all = [CLUST_all; CLUST];
            clu_str_all = [clu_str_all; clu_str];
            CluFile2load_all = [CluFile2load_all; CluFile2load] ;
            SessPath_all = [SessPath_all; SessPath];
            SessID_all = [SessID_all; SessID];
            OptoPost_all = [OptoPost_all; OptoPost];
           
            
            
        end;
    end
end
%%
ncell = ncell_all;
MouseID= MouseID_all;
Day=Day_all;
DepthZ= DepthZ_all;
ChanID = ChanID_all;
CLUST= CLUST_all;
clu_str =clu_str_all; 
CluFile2load = CluFile2load_all;
SessPath = SessPath_all;
SessID = SessID_all;
Opto_post_sess = OptoPost_all; 

CellID = [MouseID Day ChanID clu_str]; 


for ww=1:size(CellID,1) 
RowNames{ww} = CellID(ww,:)
end
%%
listcell=[]; 
listcell = table(ncell, MouseID,  Day,  ChanID, CLUST, Opto_post_sess, CluFile2load, SessPath, SessID, 'RowNames', RowNames)
save([mypath '\listcell3.mat'],'listcell');  

Tcombo =listcell(:,1:5)
save([mypath '\Tcombo.mat'],'Tcombo'); 
listcell(1:5,:)
cd(mypath)
