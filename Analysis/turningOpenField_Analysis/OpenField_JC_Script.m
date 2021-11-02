%% OpenField_JC_Script 
% First use AVIDEMUX to create the movie.mp4 : 
%       VideoOUtput   = MpegAVC(x264)
%       AudioOutput   = AAC(lav)
%       Output format = MP4v2 Muxer
%       Transform ==> CROP
%       Filter ==> Color ==> Contrast and Brightness 
%       SAVE as ==>  JCVGAT07_PRE_avid.mp4'
% JC Oct 2017 in Jaegerlab

%% 1-create MOvie FrAME (1/10) 
% PRE
movieFileName = 'D:\JC_Video\JCVGAT07_PRE_avid.mp4'
writeMovieFrames(movieFileName)

% STIM
movieFileName = 'D:\JC_Video\JCVGAT07_STIM_avid.mp4'
writeMovieFrames(movieFileName)

% POST 
movieFileName = 'D:\JC_Video\JCVGAT07_POST_avid.mp4'
writeMovieFrames(movieFileName)

%% 2- Manually Create the f0 file in each Analyzed Folder
disp('NOW CREATE THE f0.jpg MANUALLY by substracting the mouse from the frame') 
disp('copy and past backround')


%% 3- Run the Analysis (count CW and CCW rotation in 4min) 

bgFrameNo = 0 
startFrameNo = 2 
numMinutes=   2.5 %min
threshold = 0.6  % suggested between 4 and 8 (the higher the faster)
cropX=[]
cropY=[]
makeVideo=1

% PRE Manually set the 100mm distances between points.
threshold = 0.6 
movieFileName = 'D:\JC_Video\JCVGAT07_PRE_avid.mp4'
trackMouseLocomotion_JC
%%
% STIM Manually set the 100mm distances between points.
threshold = 0.75
movieFileName = 'D:\JC_Video\JCVGAT07_STIM_avid.mp4'
trackMouseLocomotion_JC

% POST Manually set the 100mm distances between points.
threshold = 0.6 
movieFileName = 'D:\JC_Video\JCVGAT07_POST_avid.mp4'
trackMouseLocomotion_JC

