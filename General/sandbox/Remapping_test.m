% Remapping_correction_Script_JC
% by Julien Catanese 9/27/2018
% 
% cd ('D:\JC_Analysis\JCVGAT07\vgat07_w25d5_z4000_25cells_200trials_tasknoopto_CAED_VM_180301_112459_novid')
% load info
% old_list = [];
% for i=1:32;
%     disp([ info.info_amp_ch(i).native_channel_name '  ' info.info_amp_ch(i).custom_channel_name])
%     old_list = [old_list; {info.info_amp_ch(i).native_channel_name}, {info.info_amp_ch(i).custom_channel_name}]
% end
% save('D:\JC_Analysis\remapping_old_list.mat', 'old_list', 'info')
% 
% cd ('D:\JC_Analysis\JCVGAT14\vgat14_w14d7_z4350_20cells_100trials_pretasknoopto_CC4F_VM_180714_202419_novid_noopto_driftperiods')
% load info
% New_list = [];
% for i=1:32;
%     disp([ info.info_amp_ch(i).native_channel_name '  ' info.info_amp_ch(i).custom_channel_name])
%     New_list = [New_list; {info.info_amp_ch(i).native_channel_name}, {info.info_amp_ch(i).custom_channel_name}]
% end
% New_list{1,2}='S3Ch8'
% New_list{7,2}='S4Ch8'
% save('D:\JC_Analysis\remapping_New_list.mat', 'New_list', 'info')
% 






%%
cd ('D:\JC_Analysis\JCVGAT07\vgat07_w25d5_z4000_25cells_200trials_tasknoopto_CAED_VM_180301_112459_novid')
load('D:\JC_Analysis\remapping_New_list.mat', 'New_list', 'info'); clear info; 
load ('info.mat');
A=dir('S*Ch*_raw.mat');
Nchanel=max(size(A));
info.info_amp_ch(1).custom_channel_name_OLDmapping=[]; 
for nch=1:Nchanel;
    disp(['current Mapping name=' info.info_amp_ch(nch).custom_channel_name]);
    B=dir([info.info_amp_ch(nch).custom_channel_name '_raw.mat']);
    clear data sr unit
    load(B.name, 'data', 'sr', 'unit');
    if isempty(info.info_amp_ch(nch).custom_channel_name_OLDmapping)
        info.info_amp_ch(nch).custom_channel_name_OLDmapping = info.info_amp_ch(nch).custom_channel_name;
    else
       disp(['old Mapping name=' info.info_amp_ch(nch).custom_channel_name_OLDmapping]);
    end
    
    for ii=1:max(size(New_list));
        if New_list{ii,1}== info.info_amp_ch(nch).native_channel_name;
            disp([New_list{ii,1} '==' info.info_amp_ch(nch).native_channel_name]);
            disp([info.info_amp_ch(nch).custom_channel_name ' will be change into ' New_list{ii,2} ]);
            
            info.info_amp_ch(nch).custom_channel_name = New_list{ii,2};
            disp(['new Mapping name=' info.info_amp_ch(nch).custom_channel_name]);
            delete(['*' info.info_amp_ch(nch).custom_channel_name '*raw*'])
            save([info.info_amp_ch(nch).custom_channel_name '_raw.mat'], 'data', 'sr', 'unit');
            save('info.mat', 'info');
            disp(['saving ' info.info_amp_ch(nch).custom_channel_name '_raw.mat  and   info.mat' ])
        end
    end
    
end


%% Rename files
% Get all text files in the current folder
% files = dir('*.lab');
% % Loop through each file
% for id = 1:1% length(files)
%     % Get the file name
%     [~, f,ext] = fileparts(files(id).name)
%     find(f,'S1Ch1')
%
% %     rename = strcat(f,'_',ext)
% %     movefile(files(id).name, rename);
% end
%

