% RenameFiles_JCscript
% to obtain M#####_yyyymmdd_sess01_rec01_SWE
% by Julien Catanese 06/09/2021

Folder = 'D:\data\Roman Exp Wheel'
cd(Folder);
file=dir([ Folder '\*_21_*.csv' ]);

for ii=1:size(file,1);
        display([file(ii).folder '\' file(ii).name]);
%         movefile([file(ii).folder '\' file(ii).name ], [file(ii).folder '\' file(ii).name(1:18)  'rec.csv']);
%         display([file(ii).folder '\' file(ii).name(1:18)  'rec.csv']);
end


