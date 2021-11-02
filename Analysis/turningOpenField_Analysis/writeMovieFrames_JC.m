function writeMovieFrames(movieFileName)

%Start the timer
timerVal = tic;

% Write out each frame to a file
% http://www.mathworks.com/help/matlab/examples/convert-between-image-sequences-and-video.html
movieObj = VideoReader(movieFileName);
frameRate = ceil(movieObj.FrameRate);
slashes = strfind (movieFileName, '\');
dots = strfind (movieFileName, '.');
outDir = movieFileName(1:slashes(end));
tmpDir = [movieFileName(1:dots(end)-1) '_analyzed'];

nFrames = movieObj.NumberOfFrames - 10; %just to be safe. to avoid EOF

%See if any missing files
if exist(tmpDir,'dir') == 0
    mkdir(tmpDir);
    disp(['Making tmp dir: ' tmpDir]);
end

imageNames = dir(fullfile(tmpDir,'f*.jpg'));
if length(imageNames) == nFrames+2
    disp('Correct number of files found. Proceeding to video analysis.');
else
    disp('Temp files missing. Regenerating tmp jpgs');
    estTimeInMin = floor((nFrames * 0.02) / 60);
    disp(['Writing ' num2str(nFrames) ' frames to disk. Estimated time to complete is ' num2str(estTimeInMin) ' minutes.']);
       
    delete([tmpDir '/*.jpg']); % remove existing files

    
    dispInterval = ceil(nFrames/100);
    jj=0;
    for ii = 1:10:nFrames
        jj = jj+1;
        %         if jj<10; prename= '0000';
        %         elseif jj>9   & jj<100  ; prename= '000';
        %         elseif jj>99  & jj<1000 ; prename= '00';
        %         elseif jj>999 & jj<10000; prename= '0';
        %         end
        
        if mod(ii, dispInterval) == 0
            disp(['   Now on frame ' num2str(ii) ' of ' num2str(nFrames)]);
        end
        %         imwrite(read(movieObj,ii),fullfile(tmpDir,sprintf(['f' prename '%d.jpg'],jj)));
        imwrite(read(movieObj,ii),fullfile(tmpDir,sprintf(['f%d.jpg'],jj)));
        
    end
end

analysisTime = toc(timerVal);
disp(['Analysis completed in ' num2str(floor(analysisTime / 60)) ' min ' num2str(mod(analysisTime,60)) ' s.' ]);
disp(['Avg time per frame was ' num2str(analysisTime / nFrames) ' s']);

clear movieObj

end
