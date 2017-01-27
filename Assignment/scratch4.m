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

imgs = zeros(size(ims,1),size(ims,2),size(ims,4));
imgsf = zeros(size(ims,1),size(ims,2),size(ims,4));
imgsf_norm = conv2(ones(size(ims,1),size(ims,2)),1/36*ones(6),'same');

for ii = 1:size(ims,4)
    imgs(:,:,ii) = rgb2gray(squeeze(ims(:,:,:,ii)));
    imgsf(:,:,ii) = conv2(double(squeeze(imgs(:,:,ii))),1/36*ones(6),'same')./imgsf_norm;
end

x_indf = min(max(abs(diff(imgsf,[],1)),[],1),[],3);
y_indf = min(max(abs(diff(imgsf,[],2)),[],2),[],3);
[~,x_ind] = max(x_indf(128:size(ims,2) - 128));
[~,y_ind] = max(y_indf(128:size(ims,1) - 128));

x_ind = x_ind + 128;
y_ind = y_ind + 128;

xshift(1) = 0;
yshift(1) = 0;
xn = conv(ones(1,size(ims,2)),ones(1,size(ims,2)));
yn = conv(ones(size(ims,1),1),ones(size(ims,1),1));
xn(1:512) = Inf;
xn(end-511:end) = Inf;
yn(1:512) = Inf;
yn(end-511:end) = Inf;

for ii = 2:size(ims,4)
    xa = conv(imgs(y_ind,:,ii),fliplr(squeeze(imgs(y_ind,:,ii-1))))./xn;
    [~,xs] = max(xa);
    xshift(ii) = xs - size(ims,2);
    vii(ii,:,1) = vii(ii-1,:,1) + xshift(ii);
    ya = conv(imgs(:,x_ind,ii),flipud(squeeze(imgs(:,x_ind,ii-1))))./yn;
    [~,ys] = max(ya);
    yshift(ii) = ys - size(ims,1);
    vii(ii,:,2) = vii(ii-1,:,2) + yshift(ii);
end



doplots = 1;
if doplots
    figure(1);
    h1 = imagesc(squeeze(ims(:,:,:,1)));
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
