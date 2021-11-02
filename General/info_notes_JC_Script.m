% info_notes_JC_Script 
% Add a info_notes field in info.mat file containing MouseID, Day and Depth.   
% written by Julien Catanese 10/2/2018


%% Add a info_notes field in info.mat file containing MouseID, Day and Depth.  
load('info.mat','info');
info.info_notes = []; clear notes;
info.info_notes.MouseID = FolderID(nf).name(1:6);
info.info_notes.Day = FolderID(nf).name(8:12);
info.info_notes.Depth = FolderID(nf).name(14:18);
info.info_notes.task = FolderID(nf).name(34:37); 

disp([info.info_notes.MouseID ' ' info.info_notes.Day ' ' info.info_notes.task ' ' info.info_notes.Depth  ]) 

%% save info
save('info.mat','info');
disp(['saving info.mat' ]);
