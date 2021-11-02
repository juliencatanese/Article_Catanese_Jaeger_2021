function  Data_allrhd  =  Get_Data_Allrhd_JCfun  (Folder)

% function Data_allrhd = Get_Data_Allrhd_JCfun(Folder)
%
% call the function: read_Intan_RHD2000_JCfun.m 
% concatenate all 1min.rhd files in Folder and output them as a struct var = 'Data_allrhd'   
%
% Catanese J. Oct 2017 JaegerLab

list_rhd = dir([Folder '\*.rhd']);
Nrhd = length(list_rhd);

for i = 1:Nrhd
    Folder = [Folder '\']; 
    File = list_rhd(i).name;
    Data_allrhd(i) = read_Intan_RHD2000_file_JCfun(Folder,File); 
end



