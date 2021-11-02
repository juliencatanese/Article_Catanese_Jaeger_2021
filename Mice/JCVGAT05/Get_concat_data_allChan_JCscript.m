% Script Get_chan_concat_data_JCscript
% 
% INPUT arg = dayFolder = 'D:\JC_Data\Acute\Awake\JC-VGAT05\day2\'
% OUTPUT = save AllChan.mat in each folder (e.g. \4000stim1\)
% Chan_data_mat = Matrice contening 32 channels data (Volt).  
% Chan_ID = Name ordered of the 32 channels in the data matrice.   
%
% Catanese J. Oct 2017 JaegerLab
clear all, close all;

dayFolder = 'D:\JC_Data\Acute\Awake\JC-VGAT05\day2\'
P = dir(dayFolder)

List_Folder={};
for nN = 1:max(size(P))
    if P(nN).isdir == 1
        List_Folder = [List_Folder P(nN).name];
    end
end
List_Folder = List_Folder(3:end)

%% LOOP ACROSS FOLDER
for nF = 1:max(size(List_Folder))
    Folder= [ P(1).folder '\' List_Folder{nF}]

    
    %% Concatene the raw data from all rhd files into a structure:
    
    Data_allrhd = Get_Data_allrhd_JC(Folder);
    
    %% Loop over 32 channels and create a Matrice of data for each channel (last column is time)
    %% Each column is a channel, [chan1 chan2, chan3, digIn1, digIn2, time]
    %% channels name are contain in a cell array that named each Matrice column.
    ChID = Data_allrhd.amplifier_channels;
    maxchan=size(ChID,2)
    
    Chan_data_mat = [];
    Chan_ID = char;
    
    for nchan=1:maxchan % usually maxchan = 32
        
        
        %% Generate concatenated data combining all files in the selected folder
        data_volt = [];
        data_digIn1 = [];
        data_digIn2 = [];
        data_time = [];
        
        for Nfile = 1:length(Data_allrhd)
            sizedata = length(Data_allrhd(Nfile).amplifier_data(1,:));
            data_volt(end+1:end+sizedata,1) = Data_allrhd(Nfile).amplifier_data(nchan,:);
        end
        
        Chan_data_mat(:,nchan)= data_volt;
        Chan_ID(:,nchan)= ChID(nchan).custom_channel_name;
    end
    
    %% SORT by SHANK and CHANNEL (S1Ch1, S1Ch2, S1Ch3, ... S2CH1, ... S4Ch8)
    % Sort by Shank only
    I1=[]; I2=[];  
    Y1=[]; Y2=[]; 
        
    sort1_Chan_data_mat = [];
    sort2_Chan_data_mat = [];
    
    sort1_Chan_ID=[];
    sort2_Chan_ID=[];
    
    ChID_I1=char; 
    ChID_I2=char; 
    ChanelID_Order = char;
    ShankID_Order = char; 
    
    [Y1,I1] = sort(Chan_ID(2,:))
    for nch = 1: maxchan
        sort1_Chan_data_mat(:,nch)=Chan_data_mat(:,I1(nch));
        sort1_Chan_ID(nch)=Chan_ID(5,I1(nch));
    end
    ChID_I1= Chan_ID(5,I1);
    ShankID_Order=Chan_ID(2,I1);
    %%
    % Sort By Channels into Shank
    sort2_Chan_data_mat=[];
    A = (Y1*100) +  sort1_Chan_ID
    [Y2, I2] = sort(A)
    for nch = 1: maxchan
        sort2_Chan_data_mat(:,nch)= sort1_Chan_data_mat(:,I2(nch));
        sort2_Chan_ID(nch)=sort1_Chan_ID(I2(nch)); % just for visual check but useless after.
    end
    ChID_I2= ChID_I1(I2);
    ChanelID_Order = ChID_I2
    
    Chan_ID(2,:)=ShankID_Order;
    Chan_ID(5,:)=ChanelID_Order;
    Chan_ID
    
    %% ADD DIGITAL INPUT1 and 2 and TIMES Chanels (Opto Stim = digI2)
    data_digIn2=[]; data_digIn1=[]; data_time = [];
    
    for Nfile = 1:length(Data_allrhd)
        sizedata = length(Data_allrhd(Nfile).amplifier_data(1,:));
        
        data_digIn1(end+1:end+sizedata,1) = Data_allrhd(Nfile).board_dig_in_data(1,:);
        data_digIn2(end+1:end+sizedata,1) = Data_allrhd(Nfile).board_dig_in_data(2,:);
        data_time(end+1:end+sizedata,1) = Data_allrhd(Nfile).t_amplifier;
    end
    
    %%
    Chan_data_mat(:,maxchan+1)= data_digIn1;
    Chan_ID(:,maxchan+1)= 'digI1';
    
    Chan_data_mat(:,maxchan+2)= data_digIn2;
    Chan_ID(:,maxchan+2)= 'digI2';
    
    Chan_data_mat(:,maxchan+3)= data_time;
    Chan_ID(:,maxchan+3)= 'times';
    
    %% SAVING
    save([Folder '\AllChan.mat'], 'Chan_data_mat', 'Chan_ID')
end
