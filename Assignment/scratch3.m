fn = fullfile('/Users/charles/Documents/MATLAB/personal/Assignment','video1.mp4');


%
coords = [17,57,558,303];
%
if ~exist('ims','var')
    v = VideoReader(fn);
    ims = read(v,[1 Inf]);
end
% % imagesc(im);aet;
% rectangle('Position',coords);

xs = [coords(1), coords(1), coords(1) + coords(3), coords(1) + coords(3), coords(1)];
ys = [coords(2), coords(2) + coords(4), coords(2) + coords(4), coords(2), coords(2)];

vsize = size(ims);
vii = zeros(vsize(4),5,2);
vii(1,:,:) = [xs.',ys.'];
vii1 = zeros(vsize(4),5,2);
vii1(1,:,:) = [xs.',ys.'];

im = squeeze(ims(:,:,:,1));
imgs = rgb2gray(im);

doplots = 0;

% Initialize features for I(1)
points = detectSURFFeatures(imgs,'MetricThreshold',64);
% points = detectBRISKFeatures(imgs,'MinContrast',.05);

[features, points] = extractFeatures(imgs, points);

features1 = features;
points1 = points;
pointsPrevious = points;
featuresPrevious = features;
pointsPreviousPrevious = pointsPrevious;
featuresPreviousPrevious = featuresPrevious;

% Initialize all the transforms to the identity matrix. Note that the
% projective transform is used here because the building images are fairly
% close to the camera. Had the scene been captured from a further distance,
% an affine transform would suffice.
tforms(size(ims,4)) = projective2d(eye(3));
tforms1(size(ims,4)) = projective2d(eye(3));

if doplots
    figure(1);
    h1 = imagesc(im);
    h2 = line(squeeze(vii(ii,1:2,1)),squeeze(vii(ii,1:2,2)));
    h3 = line(squeeze(vii(ii,2:3,1)),squeeze(vii(ii,2:3,2)));
    h4 = line(squeeze(vii(ii,3:4,1)),squeeze(vii(ii,3:4,2)));
    h5 = line(squeeze(vii(ii,4:5,1)),squeeze(vii(ii,4:5,2)));
    
end

for ii = 2:size(ims,4)
    im = squeeze(ims(:,:,:,ii));
    imgs = rgb2gray(im);
    H = hough(edge(imgs,'canny'));
    H = round(1000*H./max(H(:)));
    P = houghpeaks(H,8,'Threshold',250);
    %
    %     for ii = 1:7
    %         for ii+1:8
    %             p(end+1,:) =
    
    % Store points and features for I(n-1).
    pointsPreviousPrevious = pointsPrevious;
    featuresPreviousPrevious = featuresPrevious;
    % Store points and features for I(n-1).
    pointsPrevious = points;
    featuresPrevious = features;
    
    points = detectSURFFeatures(imgs,'MetricThreshold',64);
    %     points = detectBRISKFeatures(imgs,'MinContrast',.05);
    [features, points] = extractFeatures(imgs, points);
    
    % Find correspondences between I(n) and I(n-1).
    
    indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true);
    matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);
    indexPairs1 = matchFeatures(features, features1, 'Unique', true);
    matchedPointsPrev1 = points1(indexPairs1(:,2), :);
    
    matchedPoints = points(indexPairs(:,1), :);
    matchedPoints1 = points(indexPairs1(:,1), :);
    
    if doplots
    figure(2);
    imagesc(cat(2,squeeze(ims(:,:,:,ii-1)),squeeze(ims(:,:,:,ii))));
    aet;
    hold on;
    scatter(pointsPrevious.Location(:,1),...
        pointsPrevious.Location(:,2),'s');
    scatter(size(im,2)+points.Location(:,1),...
        points.Location(:,2),'s');
    
    plot([pointsPrevious.Location(indexPairs(:,2),1)  size(im,2)+points.Location(indexPairs(:,1),1)].',...
        [pointsPrevious.Location(indexPairs(:,2),2)  points.Location(indexPairs(:,1),2)].');
    hold off;
    end
    
    %     mp(:). = matchedPoints(:).Location
    % Estimate the transformation between I(n) and I(n-1).
    tforms(ii) = estimateGeometricTransform(matchedPoints, matchedPointsPrev,...
        'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);
    tforms1(ii) = estimateGeometricTransform(matchedPoints1, matchedPointsPrev1,...
        'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);
    
    % Compute T(1) * ... * T(n-1) * T(n)
        tforms(ii).T = tforms(ii-1).T * tforms(ii).T;
    
    vii(ii,:,:) = tforms(ii).transformPointsInverse(squeeze(vii(1,:,:)));
    vii1(ii,:,:) = tforms1(ii).transformPointsInverse(squeeze(vii1(1,:,:)));

    %     vii(ii,:,2) = vii(1,:,2);
    if doplots
        set(h1,'CData',im);
        set(h2,'XData',squeeze(vii(ii,1:2,1)),'YData',squeeze(vii(ii,1:2,2)));
        set(h3,'XData',squeeze(vii(ii,2:3,1)),'YData',squeeze(vii(ii,2:3,2)));
        set(h4,'XData',squeeze(vii(ii,3:4,1)),'YData',squeeze(vii(ii,3:4,2)));
        set(h5,'XData',squeeze(vii(ii,4:5,1)),'YData',squeeze(vii(ii,4:5,2)));
        drawnow;
    end
end

vii_orig = vii;

vii = vii_orig;
% vii = (.8*vii + .2*vii1);

% vii = vii1;

N = 20;
ff = .5*(1-cos((2*pi*(0:19))/(N-1)));
ff = ff./sum(ff(:));

normvii = conv(ones(size(vii(:,1,1))),ff,'same');
for ii = 1:5
    for jj = 1:2
        vii(:,ii,jj) = conv(vii(:,ii,jj),ff,'same')./normvii;
    end
end

doplots = 1;
if doplots
    figure(1);
    h1 = imagesc(im);
    h2 = line(squeeze(vii(ii,1:2,1)),squeeze(vii(ii,1:2,2)));
    h3 = line(squeeze(vii(ii,2:3,1)),squeeze(vii(ii,2:3,2)));
    h4 = line(squeeze(vii(ii,3:4,1)),squeeze(vii(ii,3:4,2)));
    h5 = line(squeeze(vii(ii,4:5,1)),squeeze(vii(ii,4:5,2)));
end

for ii = 2:size(ims,4)
    im = squeeze(ims(:,:,:,ii));
    imgs = rgb2gray(im);
    
    %     vii(ii,:,2) = vii(1,:,2);
    if doplots
        set(h1,'CData',im);
        set(h2,'XData',squeeze(vii(ii,1:2,1)),'YData',squeeze(vii(ii,1:2,2)));
        set(h3,'XData',squeeze(vii(ii,2:3,1)),'YData',squeeze(vii(ii,2:3,2)));
        set(h4,'XData',squeeze(vii(ii,3:4,1)),'YData',squeeze(vii(ii,3:4,2)));
        set(h5,'XData',squeeze(vii(ii,4:5,1)),'YData',squeeze(vii(ii,4:5,2)));
        drawnow;
    end
end
