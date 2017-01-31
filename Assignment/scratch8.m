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

doplots = 1;

% Initialize features for I(1)
% points = detectSURFFeatures(imgs,'MetricThreshold',512);
%     points = detectMSERFeatures(imgs,'ThresholdDelta',.5);
% points = detectBRISKFeatures(imgs,'MinContrast',.01,'MinQuality',.01);
points = detectHarrisFeatures(imgs,'FilterSize',9);

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
    h2 = line(squeeze(vii(1,1:2,1)),squeeze(vii(1,1:2,2)));
    h3 = line(squeeze(vii(1,2:3,1)),squeeze(vii(1,2:3,2)));
    h4 = line(squeeze(vii(1,3:4,1)),squeeze(vii(1,3:4,2)));
    h5 = line(squeeze(vii(1,4:5,1)),squeeze(vii(1,4:5,2)));
    
end

for ii = 2:size(ims,4)
    im = squeeze(ims(:,:,:,ii));
    imgs = rgb2gray(im)<64;
    
    % Store points and features for I(n-1).
    pointsPreviousPrevious = pointsPrevious;
    featuresPreviousPrevious = featuresPrevious;
    % Store points and features for I(n-1).
    pointsPrevious = points;
    featuresPrevious = features;
    
%     points = detectSURFFeatures(imgs,'MetricThreshold',512,'NumOctaves',4);
%     points = detectMSERFeatures(imgs,'ThresholdDelta',.5);
%         points = detectBRISKFeatures(imgs,'MinContrast',.01,'MinQuality',.01);
    points = detectHarrisFeatures(imgs,'FilterSize',27,'MinQuality',.01);
    [features, points] = extractFeatures(imgs, points);
    
    % Find correspondences between I(n) and I(n-1).
    
    indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true,'MatchThreshold',25);
    matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);
    
    matchedPoints = points(indexPairs(:,1), :);
    
    rinds = sqrt(sum((matchedPoints.Location - matchedPointsPrev.Location).^2,2))>=16;
    
    matchedPoints(rinds) = [];
    matchedPointsPrev(rinds) = [];
    
    if(length(matchedPoints)<16)
        tforms(ii).T = eye(3);
    else
%     rinds = matchedPoints.Location(:,1)>=800;
    
    if doplots
    figure(2);
    imagesc(squeeze(ims(:,:,:,ii)));
    aet;
    hold on;
    scatter(pointsPrevious.Location(:,1),...
        pointsPrevious.Location(:,2),'sg');
    scatter(points.Location(:,1),...
        points.Location(:,2),'sb');
    
    plot([matchedPointsPrev.Location(:,1)  matchedPoints.Location(:,1)].',...
        [matchedPointsPrev.Location(:,2)  matchedPoints.Location(:,2)].','s-r');
    hold off;
    end
    
    %     mp(:). = matchedPoints(:).Location
    % Estimate the transformation between I(n) and I(n-1).
    tforms(ii) = estimateGeometricTransform(matchedPoints, matchedPointsPrev,...
        'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);
    end
    Tf(ii,:,:) = tforms(ii).T;
    if(ii==72); keyboard; end;
    % Compute T(1) * ... * T(n-1) * T(n)
        tforms(ii).T = tforms(ii-1).T * tforms(ii).T;
    
    vii(ii,:,:) = tforms(ii).transformPointsInverse(squeeze(vii(1,:,:)));

    %     vii(ii,:,2) = vii(1,:,2);
    if 0
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

% N = 20;
% ff = .5*(1-cos((2*pi*(0:19))/(N-1)));
% ff = ff./sum(ff(:));
% 
% normvii = conv(ones(size(vii(:,1,1))),ff,'same');
% for ii = 1:5
%     for jj = 1:2
%         vii(:,ii,jj) = conv(vii(:,ii,jj),ff,'same')./normvii;
%     end
% end

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
