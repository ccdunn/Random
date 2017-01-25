im = rgbread(imdir('tennis.jpg'));
addpath(genpath(fullfile(homedir,'libs','edison_matlab_interface')));

seg1 = edison_wrapper(im,@RGB2Luv,'MinimumRegionArea',512);

[x,y] = ismember(squeeze(seg1(:,:,1)),unique(colvec(seg1)));
    cols = rand(max(colvec(y)),1);
    
    imseg = cols(y(:),:);
    imseg = reshape(imseg,size(im,1),size(im,2));
imagesc(imseg);
aet;

vedge = diff(imseg,[],1)~=0;
hedge = diff(imseg,[],2)~=0;

vedge = cat(1,vedge,zeros(1,size(im,2))) | cat(1,zeros(1,size(im,2)),vedge);
hedge = cat(2,hedge,zeros(size(im,1),1)) | cat(2,zeros(size(im,1),1),hedge);

col1 = (vedge | hedge);
se = strel('disk',2);
col1 = ~imdilate(col1,se);


im = rgbread(imdir('tennis.jpg'));
addpath(genpath(fullfile(homedir,'libs','edison_matlab_interface')));

seg1 = edison_wrapper(im,@RGB2Luv,'MinimumRegionArea',128);

[x,y] = ismember(squeeze(seg1(:,:,1)),unique(colvec(seg1)));
    cols = rand(max(colvec(y)),1);
    
    imseg = cols(y(:),:);
    imseg = reshape(imseg,size(im,1),size(im,2));
imagesc(imseg);
aet;

vedge = diff(imseg,[],1)~=0;
hedge = diff(imseg,[],2)~=0;

vedge = cat(1,vedge,zeros(1,size(im,2))) | cat(1,zeros(1,size(im,2)),vedge);
hedge = cat(2,hedge,zeros(size(im,1),1)) | cat(2,zeros(size(im,1),1),hedge);


col2 = ~(vedge | hedge);


mack = col1 | col2;

imagesc(mack);aet;
colormap('gray');