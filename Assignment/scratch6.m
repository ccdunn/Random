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
viiPrevPrev = zeros(vsize(4),5,2);
viiPrevPrev(1,:,:) = [xs.',ys.'];

im = squeeze(ims(:,:,:,1));
imgs = rgb2gray(im);

doplots = 0;


% Initialize all the transforms to the identity matrix. Note that the
% projective transform is used here because the building images are fairly
% close to the camera. Had the scene been captured from a further distance,
% an affine transform would suffice.

if doplots
    figure(1);
    h1 = imagesc(im);
    h2 = line(squeeze(vii(ii,1:2,1)),squeeze(vii(ii,1:2,2)));
    h3 = line(squeeze(vii(ii,2:3,1)),squeeze(vii(ii,2:3,2)));
    h4 = line(squeeze(vii(ii,3:4,1)),squeeze(vii(ii,3:4,2)));
    h5 = line(squeeze(vii(ii,4:5,1)),squeeze(vii(ii,4:5,2)));
    
end

T = zeros(size(ims,4),3,3);
T(1,:,:) = eye(3);
for ii = 1:size(ims,4)
    im = squeeze(ims(:,:,:,ii));
    imgs = rgb2gray(im);
    %
    
    if ii==2
        % Store points and features for I(n-1).
        pointsPrevious = points;
        featuresPrevious = features;
    end
    
    points = detectSURFFeatures(imgs,'MetricThreshold',64);
    %     points = detectBRISKFeatures(imgs,'MinContrast',.05);
    [features, points] = extractFeatures(imgs, points);
    
    % Find correspondences between I(n) and I(n-1).
    if ii>1
        indexPairs = matchFeatures(features, featuresPrevious, 'Unique', false);
        curmatchedPoints = points(indexPairs(:,1), :);
        matchedPoints = pointsPrevious(indexPairs(:,2), :);
        tforms(ii) = estimateGeometricTransform(curmatchedPoints, matchedPoints,...
            'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);
    else
        tforms(size(ims,4)) = projective2d(eye(3));
    end
    % Compute T(1) * ... * T(n-1) * T(n)
    T(ii,:,:) = tforms(ii).T;
    
end

T_orig = T;

T = T_orig;
N = 5;
% T = cat(1,repmat(T(1,:,:),N,1,1),T,repmat(T(1,:,:),N,1,1));

ff = .5*(1-cos((2*pi*(0:N-1))/(N-1)));
ff = ff./sum(ff(:));
ffn = conv(ones(size(T,1),1),ff,'same');

for ii = 1:3
    for jj = 1:3
        T(:,ii,jj) = conv(T(:,ii,jj),ff,'same')./ffn;
    end
end
% T = T(N+1:end-N,:,:);

Ts = repmat(max(abs(cat(1,T,T_orig))),size(T,1),1,1);
plot(T(:,:)./Ts(:,:));
hold on;
plot(T_orig(:,:)./Ts(:,:),':');
hold off;


for ii = 2:size(ims,4)
%     tforms(ii).T = squeeze(T(ii,:,:))*tforms(ii-1).T;
    tforms(ii).T = squeeze(T(ii,:,:));
    vii(ii,:,:) = tforms(ii).transformPointsInverse(squeeze(vii(1,:,:)));
end

vii_orig = vii;
vii = vii_orig;

N = 10;
ff = .5*(1-cos((2*pi*(0:19))/(N-1)));
ff = ff./sum(ff(:));
ffn = conv(ones(size(vii,1),1),ff,'same');

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
