% 
vname = 'video1';
% %coordinates of rectangle on first frame [x y width height]
coords = [17,57,558,303];
% 
% vname = 'video2';
% %coordinates of rectangle on first frame [x y width height]
% coords = [3,160,600,260];

vext = '.mp4';
vpath = fullfile(userpath,'personal/Assignment/submit',[vname vext]);

tic
videoTrack(vpath,coords);
toc