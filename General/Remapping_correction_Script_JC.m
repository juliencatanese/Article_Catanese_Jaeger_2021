% Remapping_correction_Script_JC
% by Julien Catanese 9/27/2018

%%
load('D:\JC_Analysis\remapping_New_list.mat', 'New_list', 'info'); clear info;
load ('info.mat');
A=dir('S*Ch*_raw.mat');
Nchanel=max(size(A));
info.info_amp_ch(1).custom_channel_name_OLDmapping=[];
%%
% if ~isempty(info.info_amp_ch(Nchanel).custom_channel_name_OLDmapping);
%     disp('remapping already done previously')
% else
    for nch=1:Nchanel;
        disp(['current Mapping name=' info.info_amp_ch(nch).custom_channel_name]);
        if currentFolder(25:36)== 'vgat15_w10d8' % manually corrected exceptions
            LIST= [1 5 6 8 13 15 22 23 24 26 27 29]
            nch2 = 1+LIST(nch)
            nch=nch2; 
        end
        
        B=dir([info.info_amp_ch(nch).custom_channel_name '_raw.mat']);
        clear data sr unit
        if ~isempty(B)
            disp(['loading ' B.name])
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
                    disp(['new Mapping name=' info.info_amp_ch(nch).custom_channel_name])
                    delete(['*' info.info_amp_ch(nch).custom_channel_name '*raw*'])
                    save([info.info_amp_ch(nch).custom_channel_name '_raw.mat'], 'data', 'sr', 'unit');
                    
                    
                    
                    save('info.mat', 'info');
                    disp(['saving ' info.info_amp_ch(nch).custom_channel_name '_raw.mat  and   info.mat' ])
                end
            end
        else
            disp('this ch.raw.mat file have been removed, probably had no spikes')
            info.info_amp_ch(nch).custom_channel_name_OLDmapping = info.info_amp_ch(nch).custom_channel_name;
            info.info_amp_ch(nch).custom_channel_name = 'remov'
        end
    end
% end
%% 
A=dir('S*Ch*_raw.mat');
Nchanel=max(size(A));
if Nchanel ~= max(size(info.info_amp_ch))
    disp('ATTENTION YOU NEED TO REMOVED ONE CHANNEL: ')
    for nch=1:Nchanel;
        nmatch =0; 
        for jj=1:max(size(info.info_amp_ch));
            if A(nch).name(1:5)==info.info_amp_ch(jj).custom_channel_name
                nmatch=1; 
            end
        end
        if nmatch ==0
            disp (['deleting '  A(nch).name])
            delete(['*'  A(nch).name(1:5) '*'])
        else
%             disp(['ok for ' A(nch).name])
        end
    end
else 
    disp('GOOD: Nchanel correspond to info.info_amp_ch')
end

