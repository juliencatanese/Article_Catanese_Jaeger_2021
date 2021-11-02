%% Mouse Tracking software

% Adapted from Gomez-Marin et al. (2012). "Automated Tracking of Animal
% Posture and Movement during Exploration and Sensory Orientation
% Behaviors" PloS One

% function output = trackMouseLocomotion(movieFileName, bgFrameNo, startFrameNo, numMinutes, threshold, cropX, cropY, makeVideo)
% output = trackMouseLocomotion('D:\JC_Video\JCVGAT07_PRE_avid.mp4', 0, 3, 5, 0.55, [], [], 1)
% threshold bettween 0.4 and 0.6 ==> 0.55 should work well.  
% bgFrameNo = 0 => mean that f0 is the background frame 
% startFrameNo = 3 ==> Start at frame 3 
% numMinutes = 5 ==> 5min analysis. 
% makeVideao = 1 or 0 ==> yes or no. 
% crop = [] ==> no need since we already crop in avidemux. 
% What is the distance from point A to point B (in mm)? 100 ==> 100mm between the 2 points


% cd movieFileName
% A=dir; AA=
% 
% A1; for i=1:2:size(A1);  delete(A1(i).name); end
% 
% A2=dir; for i=1:2:size(A2);  delete(A2(i).name); end
% 
% A3=dir; for i=3:size(A3);  movefile([A3(i).name], ['f' num2str(i) '.jpg']); end
% 
% cd C:\Users\JCATANE\Documents\MATLAB

% Script trackMouseLocomotion(movieFileName, bgFrameNo, startFrameNo, numMinutes, threshold, cropX, cropY, makeVideo)
%% Variables def:  
% movieFileName = 'D:\JC_Video\JCVGAT07_POST_avid.mp4'
% bgFrameNo=0 
% startFrameNo=2 
% numMinutes=2.5
% threshold=0.6 
% cropX=[]
% cropY=[]
% makeVideo=1


%Start the timer
timerVal = tic;

% Write out each frame to a file
% http://www.mathworks.com/help/matlab/examples/convert-between-image-sequences-and-video.html
movieObj = VideoReader(movieFileName);
frameRate = ceil(movieObj.FrameRate);
slashes = strfind (movieFileName, '/');
dots = strfind (movieFileName, '.');
%outDir = movieFileName(1:slashes(end));
tmpDir = [movieFileName(1:dots(end)-1) '_analyzed'];
outputVideoName = [movieFileName(1:dots(end)-1) '_analyzed.avi'];

if exist(outputVideoName,'file') > 0
    delete(outputVideoName);
end

endFrameNo = (1/10)*ceil(movieObj.FrameRate)*60*numMinutes + startFrameNo;
nFrames = endFrameNo - startFrameNo;

disp(['Analyzing movie ''' movieFileName '''...']);
delete([tmpDir '/af*.jpg']); % remove existing analysis files, if any

%See if any missing files
% if exist(tmpDir,'dir') == 0
%     mkdir(tmpDir);
%     disp(['Making tmp dir: ' tmpDir]);
% end
imageNames = dir(fullfile(tmpDir,'f*.jpg'));
if length(imageNames) >= nFrames/10+2
     disp('Tmp files found. Proceeding to video analysis.');
else
    error('Frame files missing. Please regenerate using writeMovieFrames()');
end
%     disp('Temp files missing. Regenerating tmp jpgs');
%     disp(['Writing ' num2str(nFrames) ' frames to disk. This may take a few minutes.']);
%     delete([tmpDir '/*.jpg']); % remove existing files
%     dispInterval = round(nFrames/20);
%     imwrite(read(movieObj,bgFrameNo),fullfile(tmpDir,sprintf('f%d.jpg',bgFrameNo)));
%     for ii = startFrameNo:endFrameNo
%         if mod(ii, dispInterval) == 0
%             disp(['   Now on frame ' num2str(ii) ' of ' num2str(endFrameNo)]);
%         end
%         imwrite(read(movieObj,ii),fullfile(tmpDir,sprintf('f%d.jpg',ii)));
%     end
% end



clear movieObj

% Background subtraction (a simple way to detect the animal)
ibgOrig = imread(fullfile(tmpDir,sprintf('f%d.jpg',bgFrameNo)));
%ibg = ibgOrig(cropX,cropY,3);
ibg = ibgOrig;
ibg_bw = im2bw(ibg,threshold);

%Bring up the first image and allow the user to select the scale bar marks
%figure; imshow(ibgOrig(cropX,cropY));
figure; imshow(ibg);
title('Click the two dots representing your scale bar and then enter in their measurement into the Matlab Command Window')
% In the sample data, the center of the droplet was at the South-East well
% from the middle little dark mark of the 96-well plate
set(gca,'xtick',[],'ytick',[])
[X1,Y1]=ginput(1);
[X2,Y2]=ginput(1);
close all
refDistanceInMM = input('What is the distance from point A to point B (in mm)?');

% Initialize our analyzed output
output.xctime=zeros(1,nFrames+1);
output.yctime=zeros(1,nFrames+1);
output.orientationtime=zeros(1,nFrames+1);
output.netRotations = 0;
output.totalCWRotations = 0;
output.totalCCWRotations = 0;
output.distancetime = zeros(1,nFrames+1);
output.totalDistance = 0;
output.outputVideoFileName = outputVideoName;
output.angularTotal = zeros(1,nFrames+1); 
startingOrientation = NaN;
% figure
%Analyze it!
angularTotal = 0;
dispInterval = round(nFrames/20);
disp(['Analyzing frames ' num2str(startFrameNo) ' through ' num2str(endFrameNo) ' ...']);
for i = startFrameNo:1:endFrameNo
    
    if mod(i-startFrameNo, dispInterval) == 0
        disp(['   Now processing frame ' num2str(i) '.']);
    end
    
    statIndex = i - startFrameNo + 1;
    
    % Read in the ith frame and subtract it from the background
    [i1,map] = imread(fullfile(tmpDir,sprintf('f%d.jpg',i)));
    %i1 = i1(cropX,cropY);
    i1_bw = im2bw(i1,threshold);
%     i2_bw = im2bw((i1-1).*-1,0.7);
%     i1_bw2 = im2bw(i1,0.6);
%     i1_bw3 = im2bw(i1,0.7);
%     subplot(1,2,1), imagesc(i1_bw), caxis([0,1]), axis square
%     subplot(1,2,2), imagesc(i2_bw), caxis([0,1])
%     pause(0.01)
    
    i1_bw_sub = ibg_bw - i1_bw;
    
    % Smooth it
    h=fspecial('gaussian', 3, 1);
    i1_bw_sub=imfilter(i1_bw_sub,h,'replicate');
    
    % Fill it
    i1_bw_sub=imfill(i1_bw_sub,'holes');
    
    % Clear border trick
   % i1_bw_sub=imclearborder(i1_bw_sub);
    
    % Delete small objects
    i1_bw_sub=bwareaopen(i1_bw_sub,500);
    
    % Region detection
    [L,numobj] = bwlabel(i1_bw_sub,8);
    regprops=regionprops(L,'Centroid','Area','Orientation','MajorAxisLength','MinorAxisLength'); % 'Perimeter'
    
    % Take only the big object
    bigobject1=find([regprops.Area]==max([regprops.Area]));
    bigobject=bigobject1(1);
    
    % The first and most basic coordinate: its centroid position
    centr=regprops(bigobject).Centroid;
    xc=centr(1);
    yc=centr(2);
    output.xctime(statIndex) = xc;
    output.yctime(statIndex) = yc;
    
    %Calculate distance from last center
    if i ~= startFrameNo
        output.distance(statIndex) = sqrt( (xc - output.xctime(statIndex - 1))^2 + (yc - output.yctime(statIndex - 1))^2 );
    else
        output.distance(statIndex) = 0;
    end
    
    %Determine the orientation of the biggest region
    currentOrientation = regprops(bigobject1).Orientation;
    output.orientationtime(statIndex) = currentOrientation;
    
    %Determine if the animal has rotated
    rotationLabel = '';
    if i == startFrameNo
        startingOrientation = currentOrientation;
        lastOrientation = NaN;
    else
        lastOrientation = output.orientationtime(statIndex -1);
        
        %Need to deal with the discontinuity around +/- 90 degrees
        if lastOrientation < -20 && currentOrientation > 70
            currentOrientation = currentOrientation - 180;
        elseif lastOrientation > 70 && currentOrientation < -20;
            currentOrientation = currentOrientation + 180;
        end
        
        angularTotal = angularTotal + (currentOrientation - lastOrientation);
        output.angularTotal = angularTotal;
        
        if currentOrientation > startingOrientation && startingOrientation > lastOrientation && abs(angularTotal) >= 360
            output.totalCCWRotations = output.totalCCWRotations + 1;
            angularTotal = currentOrientation - startingOrientation;
            rotationLabel = 'CCW Rotation!';
        elseif currentOrientation < startingOrientation && startingOrientation < lastOrientation && abs(angularTotal) >= 360
            output.totalCWRotations = output.totalCWRotations + 1;
            angularTotal = currentOrientation - startingOrientation;
            rotationLabel = 'CW Rotation!';
        end
    end
    
    if makeVideo
        %Write out analyzed output frame to tmp directory
        figHandle = figure('Visible','off'); hold on;
        %figHandle = figure; hold on;
        imshow(i1,map,'InitialMagnification', 100,'Border','tight');
        %image(i1);
        %colormap(gray);
        axis off;
        hold on
        
        %Mark the center of the centroid with a dot
        plot(xc,yc,'*y')
        
        %Draw lines for the major and minor axes along the orientation
        yAxis = [1 -1] * regprops(bigobject).MajorAxisLength/2 * sin(currentOrientation * (pi / 180)) + yc;
        xAxis = [-1 1] * regprops(bigobject).MajorAxisLength/2 * cos(currentOrientation * (pi / 180)) + xc;
        plot(xAxis,yAxis,'-y','LineWidth',1);
        yAxis = [1 -1] * regprops(bigobject).MinorAxisLength/2 * sin((currentOrientation+90) * (pi / 180)) + yc;
        xAxis = [-1 1] * regprops(bigobject).MinorAxisLength/2 * cos((currentOrientation+90) * (pi / 180)) + xc;
        plot(xAxis,yAxis,'--y','LineWidth',1);
        
        %Draw the perimeter
        iperim=bwperim(i1_bw_sub);
        centerIndp=find(iperim>0);
        [py,px]=ind2sub(size(iperim),centerIndp);
        contour=bwtraceboundary(iperim, [py(1) px(1)], 'N');
        plot(contour(:,2), contour(:,1), 'linewidth', 1, 'color', 'b');
        
        
        %Make a label on the bottom of the image
        newLabel = ['F#' num2str(i) ' Angle=' num2str(round(currentOrientation)) ' LastAngle= ' num2str(round(lastOrientation)) '; CW#: ' num2str(output.totalCWRotations) ' CCW#: ' num2str(output.totalCCWRotations) ' AngTot: ' num2str(round(angularTotal*10)/10) '  ' rotationLabel];
        text('units','pixels','position',[0 10],'fontsize',14,'string',newLabel,'Color','r')
        
        drawnow
        
        % And save the image for every processed frame
        print('-f1','-djpeg',strcat(tmpDir,'/af',num2str(i),'.jpg'));
        
        close(figHandle);
    end
    
end

%Determine output.distancePerMinute
output.distancePerMinute = zeros(1,floor(numMinutes));
cgTime = 1;
for i = 1:floor(numMinutes)
    output.distancePerMinute(i) = sum(output.distance(frameRate*60*(i-1)/cgTime+1 : frameRate*60*i/cgTime));
    disp(['Distance in Minute #' num2str(i) ' was ' num2str(round(output.distancePerMinute(i) * (sqrt( (X1-X2)^2 + (Y1-Y2)^2 ) / refDistanceInMM))) ' mm']);
end

output.totalDistance = sum(output.distance) * (sqrt( (X1-X2)^2 + (Y1-Y2)^2 ) / refDistanceInMM);
disp(['Total Distance in ' num2str(numMinutes) ' was ' num2str(output.totalDistance) ' mm']);

output.netRotations = output.totalCCWRotations - output.totalCWRotations;
disp(['Total CW rotations in ' num2str(numMinutes) ' min was ' num2str(output.totalCWRotations) ' rotations.']);
disp(['Total CCW rotations in ' num2str(numMinutes) ' min was ' num2str(output.totalCCWRotations) ' rotations.']);
disp(['Net rotations in ' num2str(numMinutes) ' min was ' num2str(output.netRotations) ' rotations.']);


%Make the output video and remove the temporary directory
if makeVideo
    disp(['Writing output video from analyzed frames: ''' outputVideoName '''']);
    imageNames = dir(fullfile(tmpDir,'af*.jpg'));
    imageNames = {imageNames.name}';
    imageStrings = regexp([imageNames{:}],'(\d*)','match');
    imageNumbers = str2double(imageStrings);
    [~,sortedIndices] = sort(imageNumbers);
    sortedImageNames = imageNames(sortedIndices);
    outputVideo = VideoWriter(outputVideoName);
    outputVideo.FrameRate = frameRate;
    open(outputVideo);
    for ii = 1:length(sortedImageNames)
        img = imread(fullfile(tmpDir,sortedImageNames{ii}));
        writeVideo(outputVideo,img);
    end
    close(outputVideo);
end

%How'd we do on time?
analysisTime = toc(timerVal);
disp(['Analysis completed in ' num2str(floor(analysisTime / 60)) ' min ' num2str(mod(analysisTime,60)) ' s.' ]);


