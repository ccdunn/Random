% fn = fullfile('/Users/charles/Documents/MATLAB/personal/Assignment','video1.mp4');
% v = VideoReader(fn);
% 
% coords = [17,57,558,303];
% 
% ims = read(v,[1 Inf]);
% % imagesc(im);aet;
% rectangle('Position',coords);

vsize = size(ims);

im = squeeze(ims(:,:,:,1));
imgs = rgb2gray(im);
img = imgradient(imgs);
edgim = edge(imgs,'canny');
H = hough(edgim);
tsum = sum(H,1);
xedge = sum(img,1);
yedge = sum(img,2);

im1 = im;
img21 = imgs;
img1 = img;
edgim1 = edgim;
H1 = H;
tsum1 = tsum;
xedge1 = xedge;
yedge1 = yedge;
    P1 = houghpeaks(H1,8);

tcnorm = conv(ones(1,length(tsum)),fliplr(ones(1,length(tsum))),'same');

% out = zeros([vsize(1:2) 1 vsize(4)]);

for ii = 1:size(ims,4)
    im = squeeze(ims(:,:,:,ii));
    imgs = rgb2gray(im);
    img = imgradient(imgs);
    edgim = edge(imgs,'canny');
    H = hough(edgim);
    P = houghpeaks(H,8);
    tsum = sum(H,1);
    tc = conv(tsum1,fliplr(tsum),'same')./tcnorm;
    [~,mind] = max(tc);
    toffset = mind - length(tsum)/2
    xedge = sum(img,1);
    xc = conv(xedge1,fliplr(xedge),'same');
    [~,mind] = max(xc);
    xoffset = mind - length(xedge)/2;
    yedge = sum(img,2);
    yc = conv(yedge1,fliplr(yedge),'same');
    [~,mind] = max(yc);
    yoffset = mind - length(yedge)/2;
    
%     if ii==1
%         h1 = imagesc(im);
%         h2 = rectangle('Position',coords);
%     else
%         set(h1,'CData',im);
%         set(h2,'Position',coords + [-xoffset yoffset 0 0]);
%     end
%     drawnow;
end