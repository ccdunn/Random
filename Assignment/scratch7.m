
vname = 'video1';
vext = '.mp4';
vpath = fullfile(userpath,'personal/Assignment',[vname vext]);

%coordinates of rectangle on first frame [x y width height]
coords = [17,57,558,303];

%
v = VideoReader(vpath);
ims = read(v,[1 Inf]);
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

No = 64;
Nn = 64;
Ns = 2;
Nt = .1;
[im1, im1g] = imgradient(squeeze(imgsf(No:end-No,No:end-No,1)));
im1gx = im1g(im1(:)>100);
im1gx = hist(mod(im1gx+90,180)-90,-90:Nt:90);
im1gx = im1gx./sum(im1gx);
xst = 0;
yst = 0;
for ff = 2:size(ims,4)
    [im2, im2g] = imgradient(squeeze(imgsf(No:end-No,No:end-No,ff)));
    
    Nm = 4;
    
    dneigh = nan(65,65);
    dneigh(yst+Nn/2+1,xst+Nn/2+1) = sum(colvec((...
        im1(1-min(0,yst):size(im1,1)+min(0,-yst),1-min(0,xst):size(im1,2)+min(0,-xst)) -...
        im2(1-min(0,-yst):size(im2,1)+min(0,yst),1-min(0,-xst):size(im2,2)+min(0,xst))...
        ).^2));
    dmin = nanmin(dneigh(:));
    dminprev = Inf;
    xst = min(Nn/2-Ns,max(-(Nn/2-Ns),xst));
    yst = min(Nn/2-Ns,max(-(Nn/2-Ns),yst));
    while(1)
        dminprev = nanmin(dneigh(:));
        foundmin = 0;
        for ii = yst-Ns:yst+Ns
            for jj = xst-Ns:xst+Ns
                if(ii==0 && jj==0 || (abs(ii - yst)+abs(jj - xst))>=2*Ns);
                    continue;
                end
                if(isnan(dneigh(ii+Nn/2+1,jj+Nn/2+1)))
                    dneigh(ii+Nn/2+1,jj+Nn/2+1) = sum(colvec((...
                        im1(1-min(0,ii):size(im1,1)+min(0,-ii),1-min(0,jj):size(im1,2)+min(0,-jj)) -...
                        im2(1-min(0,-ii):size(im2,1)+min(0,ii),1-min(0,-jj):size(im2,2)+min(0,jj))...
                        ).^2));
                    if(dneigh(ii+Nn/2+1,jj+Nn/2+1)<=dminprev)
                        foundmin = 1;
                    end
                end
            end
        end
        [~,mind] = nanmin(dneigh(:));
        [yst,xst] = ind2sub(size(dneigh),mind);
        xst = xst - (Nn/2+1);
        yst = yst - (Nn/2+1);
        ys(ff) = yst;
        xs(ff) = xst;
        if(~foundmin || abs(xst)>=Nn/2-Ns || abs(yst)>=Nn/2-Ns)
            break;
        end
    end
    
    %     imagesc(dneigh);aet;
    %     drawnow;
    
    vii(ff,:,1) = vii(ff-1,:,1) + xs(ff);
    vii(ff,:,2) = vii(ff-1,:,2) + ys(ff);
    
    im2gx = im2g(im2(:)>100);
    im2gx = hist(mod(im2gx+90,180)-90,-90:Nt:90);
    im2gx = im2gx./sum(im2gx);
    
    %     plot(im2gx);hold on;plot(im1gx);hold off;drawnow;
    
    as = (-15/Nt:15/Nt);
    d = zeros(1,length(as));
    for ii = 1:length(as)
        d(ii) = sum(im1gx.*circshift(im2gx,as(ii)));
    end
    [~,mind] = max(d);
    
    theta(ff) = as(mind)*Nt;
    R = [cosd(theta(ff)) -sind(theta(ff));sind(theta(ff)) cosd(theta(ff))];
    for ii = 1:5
        vii(ff,ii,1:2) = R*(squeeze(vii(ff,ii,1:2)) - [size(im,1)/2; size(im,2)/2]) + [size(im,1)/2; size(im,2)/2];
    end
    %
    im1 = im2;
    im1g = im2g;
    im1gx = im2gx;
end

doplots = 1;
if doplots
    figure(1);
    h1 = imagesc(squeeze(ims(:,:,:,1)));
    h2 = line(squeeze(vii(1,1:2,1)),squeeze(vii(1,1:2,2)));
    h3 = line(squeeze(vii(1,2:3,1)),squeeze(vii(1,2:3,2)));
    h4 = line(squeeze(vii(1,3:4,1)),squeeze(vii(1,3:4,2)));
    h5 = line(squeeze(vii(1,4:5,1)),squeeze(vii(1,4:5,2)));
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

vo = VideoWriter(fullfile(userpath,'personal/Assignment',[vname '_out']));
vo.FrameRate = v.FrameRate;  % Default 30
open(vo);

for ii = 1:size(ims,4)
    im = squeeze(ims(:,:,:,ii));
    vt = squeeze(vii(ii,:,:));
    M = [vt(1,1:2) vt(2,1:2);...
        vt(2,1:2) vt(3,1:2);...
        vt(3,1:2) vt(4,1:2);...
        vt(4,1:2) vt(5,1:2)];
    imo = insertShape(im,'Line',M,'LineWidth',5,'Color',[1 0 0]);
    writeVideo(vo, imo);
end
close(vo);


