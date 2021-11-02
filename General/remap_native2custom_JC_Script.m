% Remap_info_ampCh_custom_name_JC_Script
% to use right after dat2mat_JC_Script.m
% written by Julien Catanese 10/02/2018
% last updated JC 10/02/2018

%% Save the right Mapping in 'D:\JC_Analysis\remapping_New_list.mat'
AA = dir('D:\JC_Analysis\remapping_New_list.mat')
if isempty(AA)
    cd ('D:\JC_Analysis\JCVGAT14\vgat14_w14d7_z4350_20cells_100trials_pretasknoopto_CC4F_VM_180714_202419_novid_noopto_driftperiods')
    load info
    New_list = [];
    for i=1:32;
        disp([ info.info_amp_ch(i).native_channel_name '  ' info.info_amp_ch(i).custom_channel_name])
        New_list = [New_list; {info.info_amp_ch(i).native_channel_name}, {info.info_amp_ch(i).custom_channel_name}]
    end
    New_list{1,2}='S3Ch8'
    New_list{7,2}='S4Ch8'
    save('D:\JC_Analysis\remapping_New_list.mat', 'New_list', 'info')
end

%% Rename custom channel in info.mat file 
load('D:\JC_Analysis\remapping_New_list.mat', 'New_list');
clear info;
load('info.mat','info');
currentFolder =pwd; 
AAA = dir([pwd '\A*_raw.mat']);
Nchanel=max(size(AAA));
info.info_amp_ch(1).custom_channel_name_OLDmapping=[];

for nch=1:Nchanel;
        if currentFolder(25:36)== 'vgat15_w10d8' % manually corrected exceptions
            LIST= [1 5 6 8 13 15 22 23 24 26 27 29]
            nch2 = 1+LIST(nch)
            nch=nch2; 
        end
    disp(['Native name = '  info.info_amp_ch(nch).native_channel_name           ' current Mapping name = ' info.info_amp_ch(nch).custom_channel_name ])
    disp(['loading ' info.info_amp_ch(nch).native_channel_name '_raw.mat'])
    load([info.info_amp_ch(nch).native_channel_name '_raw.mat'])  
    disp(['loading ' info.info_amp_ch(nch).native_channel_name])
    for ii=1:max(size(New_list));
        if  info.info_amp_ch(nch).native_channel_name == New_list{ii,1} 
            disp([info.info_amp_ch(nch).native_channel_name   '==' New_list{ii,1}]); % New_list{ii,1} = native name 
            disp([info.info_amp_ch(nch).custom_channel_name ' will be change into ' New_list{ii,2} ]); % New_list{ii,2} = corresponding correct custom name 
            
            info.info_amp_ch(nch).custom_channel_name = New_list{ii,2};
            disp(['new Mapping name=' info.info_amp_ch(nch).custom_channel_name]);
            
            save([info.info_amp_ch(nch).custom_channel_name '_raw.mat'], 'data', 'sr', 'unit');
            disp(['saving ' info.info_amp_ch(nch).custom_channel_name '_raw.mat' ])
            
        end
    end
        
    if isempty(info.info_amp_ch(nch).custom_channel_name_OLDmapping)
        info.info_amp_ch(nch).custom_channel_name_OLDmapping = info.info_amp_ch(nch).custom_channel_name;
    else
        disp(['old Mapping name=' info.info_amp_ch(nch).custom_channel_name_OLDmapping]);
    end
end

%% save info
save('info.mat','info');
disp(['saving info.mat' ]);

