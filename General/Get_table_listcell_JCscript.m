% Get_table_listcell_JCcript
% create and save a table that contain all cells informations
% type: load('D:\JC_Analysis\listcell.mat','listcell')
% type:listcell
% written by Julien Catanese 12/05/2018

clear all; 

ncell=0;
ncell_all = [] ;
MouseID_all = [];
Day_all =[];
DepthZ_all = [];
ChanID_all = [];
CLUST_all = [];
clu_str_all = [];
CluFile2load_all = [] ;
SessPath_all = [];
SessID_all = [];

% cd('D:\JC_Analysis');
AnalysisFolder = 'C:\Users\catan\Documents\EMORY\JC_Analysis'
cd(AnalysisFolder)
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
    
    clust_file = dir('times_*S*Ch*_sub.mat');   ncluf=max(size(clust_file));
    for nchan = 1:ncluf %  channel loop
        ChanID=clust_file(nchan).name(7:11);
        load(clust_file(nchan).name);
        Nclust = max(cluster_class(:,1));
        for CLUST= 1:Nclust %Cluster loop within each channel      
            ncell=ncell+1; 
            
            clu_str = ['clu#' num2str(CLUST)]; 
            CluFile2load =  [clust_file(nchan).name];
            
            disp ([MouseID ' ' Day  ' ' DepthZ ' ' ChanID ' clu#' num2str(CLUST)])
          
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

CellID = [MouseID Day DepthZ ChanID clu_str]; 

for ww=1:size(CellID,1) 
RowNames{ww} = CellID(ww,:)
end

listcell = table(ncell, MouseID,  Day, DepthZ, ChanID, CLUST, CluFile2load, SessPath, SessID, 'RowNames', RowNames)
save([AnalysisFolder '\listcell3.mat'],'listcell')

listcell(1:5,:)

