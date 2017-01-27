
vname = 'video1';
vext = '.mp4';
vpath = fullfile(userpath,'personal/Assignment',[vname vext]);

%coordinates of rectangle on first frame [x y width height]
coords = [17,57,558,303];

% vname = 'video2';
% vext = '.mp4';
% vpath = fullfile(userpath,'personal/Assignment',[vname vext]);
% 
% %coordinates of rectangle on first frame [x y width height]
% coords = [3,160,600,260];


videoTrack(vpath,coords);


