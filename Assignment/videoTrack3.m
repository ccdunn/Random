% function videoTrack3(vpath,coords)


vname = 'video1';
vext = '.mp4';
vpath = fullfile(userpath,'personal/Assignment',[vname vext]);

%coordinates of rectangle on first frame [x y width height]
coords = [17,57,558,303];

%load video frames
v = VideoReader(vpath);
ims = read(v,[1 Inf]);

%input video size
Nv_h = size(ims,1);
Nv_w = size(ims,2);
Nv_c = size(ims,3);
Nv_f = size(ims,4);

%initial corners of rectangle
x_i = [coords(1); coords(1); coords(1) + coords(3); coords(1) + coords(3)];
y_i = [coords(2); coords(2) + coords(4); coords(2) + coords(4); coords(2)];

%initialize coordinates for each frame
p_rect = zeros(Nv_f,4,2);
p_rect(1,:,:) = [x_i,y_i];

%convert each frame to grayscale
imgs = zeros(Nv_h,Nv_w,Nv_f);
imgsf = zeros(Nv_h,Nv_w,Nv_f);


%smoothing helps with some noise issues
N = 8;
ff = .5*(1-cos((2*pi*(0:N-1))/(N-1)));
ff = ff./sum(ff(:));
ff = ff.'*ff;
imgsf_norm = conv2(ones(Nv_h,Nv_w),ff,'same');
for ii = 1:Nv_f
    imgs(:,:,ii) = rgb2gray(squeeze(ims(:,:,:,ii)));
    imgsf(:,:,ii) = conv2(squeeze(imgs(:,:,ii)),ff,'same')./imgsf_norm;
end  


%convert each frame to grayscale
E = zeros(Nv_h,Nv_w,Nv_f);

%smoothing helps with some noise issues
for ii = 1:Nv_f
    E(:,:,ii) = edge(squeeze(imgsf(:,:,ii)),'canny',.05);
end  
    

%convert each frame to grayscale
clear H;
P = zeros(Nv_f,2,2);

rho_res = 5;

%smoothing helps with some noise issues
for ii = 1:Nv_f
    [H(:,:,ii),Th,Rh] = hough(squeeze(E(:,:,ii)),'rhoresolution',rho_res,'Theta',-90:89);    
    P(ii,:,:) = houghpeaks(squeeze(H(:,:,ii)),2,'Threshold',eps);
end  
    

No = 64;%buffer region to avoid objects going out of frame
Nn = 64;%number of offsets in each direction to check
Ns = 2;%size of local minimum check in each direction
Nt = .1;%angle resolution in degrees for rotation search

%calculate metrics for first frame
[im1, im1g] = imgradient(squeeze(imgsf(No:end-No,No:end-No,1)));
im1gx = im1g(im1(:)>100);
im1gx = hist(mod(im1gx+90,180)-90,-90:Nt:90);
im1gx = im1gx./sum(im1gx);
xst = 0;
yst = 0;

rho_res = 1;

xs = zeros(1,Nv_f);
ys = zeros(1,Nv_f);
%iterate over frames
for ff = 2:Nv_f
    [im2, im2g] = imgradient(squeeze(imgsf(No:end-No,No:end-No,ff)));
    
    Ht = squeeze(H(:,:,ff));
    Ht(:,abs(mod(Th-Th(P(ff-1,1,2))+90,180)-90)>5) = 0;
    Ht((abs(Rh-Rh(P(ff-1,1,1)))>32),:) = 0;
    Pt = houghpeaks(Ht,1,'Threshold',eps);
    
    %start search for SSD minimum at previous shift
    ssd = nan(65,65);
    ssd(yst+Nn/2+1,xst+Nn/2+1) = sum(colvec((...
        im1(1-min(0,yst):size(im1,1)+min(0,-yst),1-min(0,xst):size(im1,2)+min(0,-xst)) -...
        im2(1-min(0,-yst):size(im2,1)+min(0,yst),1-min(0,-xst):size(im2,2)+min(0,xst))...
        ).^2));
    xst = min(Nn/2-Ns,max(-(Nn/2-Ns),xst));
    yst = min(Nn/2-Ns,max(-(Nn/2-Ns),yst));
    %gradient descent
    while(1)
        dminprev = nanmin(ssd(:));
        foundmin = 0;
        for ii = yst-Ns:yst+Ns
            for jj = xst-Ns:xst+Ns
                if((abs(ii - yst)+abs(jj - xst))>=2*Ns)
                    continue;
                end
                if(isnan(ssd(ii+Nn/2+1,jj+Nn/2+1)))
                    ssd(ii+Nn/2+1,jj+Nn/2+1) = sum(colvec((...
                        im1(1-min(0,ii):size(im1,1)+min(0,-ii),1-min(0,jj):size(im1,2)+min(0,-jj)) -...
                        im2(1-min(0,-ii):size(im2,1)+min(0,ii),1-min(0,-jj):size(im2,2)+min(0,jj))...
                        ).^2));
                    if(ssd(ii+Nn/2+1,jj+Nn/2+1)<=dminprev)
                        foundmin = 1;
                    end
                end
            end
        end
        [~,mind] = nanmin(ssd(:));
        [yst,xst] = ind2sub(size(ssd),mind);
        xst = xst - (Nn/2+1);
        yst = yst - (Nn/2+1);
        ys(ff) = yst;
        xs(ff) = xst;
        if(~foundmin || abs(xst)>=Nn/2-Ns || abs(yst)>=Nn/2-Ns)
            break;
        end
    end
    
    dnm(ff,:,:) = ssd;

%     figure(1);
%     imagesc(-32:32,-32:32,ssd);aet;
%     xlabel('Horizontal shift (pixels)');
%     ylabel('Vertical shift (pixels)');
%     title('SSD');
%     colorbar;
%     saveas(1,'SSD_plot.png');
    
    %shift rectangle coordinates based on determined frame shift
    p_rect(ff,:,1) = p_rect(ff-1,:,1) + xs(ff);
    p_rect(ff,:,2) = p_rect(ff-1,:,2) + ys(ff);
    
    
    im2gx = im2g(im2(:)>100);
    im2gx = hist(mod(im2gx+90,180)-90,-90:Nt:90);
    im2gx = im2gx./sum(im2gx);
    
    %     plot(im2gx);hold on;plot(im1gx);hold off;drawnow;
    
    as = (-5/Nt:5/Nt);
    ssd = zeros(1,length(as));
    for ii = 1:length(as)
        ssd(ii) = sum(im1gx.*circshift(im2gx,as(ii)));
    end
    [~,mind] = max(ssd);
    
    theta(ff) = Th(Pt(2)) - Th(P(ff-1,1,2));%as(mind)*Nt;
    R = [cosd(theta(ff)) -sind(theta(ff));sind(theta(ff)) cosd(theta(ff))];
    for ii = 1:4
        p_rect(ff,ii,1:2) = R*(squeeze(p_rect(ff,ii,1:2)) - [Nv_h/2; Nv_w/2]) + [Nv_h/2; Nv_w/2];
    end
    %
    im1 = im2;
    im1g = im2g;
    im1gx = im2gx;
end

% doplots = 1;
% if doplots
%     figure(1);
%     h1 = imagesc(squeeze(ims(:,:,:,1)));
%     h2 = line(squeeze(vii(1,1:2,1)),squeeze(vii(1,1:2,2)));
%     h3 = line(squeeze(vii(1,2:3,1)),squeeze(vii(1,2:3,2)));
%     h4 = line(squeeze(vii(1,3:4,1)),squeeze(vii(1,3:4,2)));
%     h5 = line(squeeze(vii(1,4:5,1)),squeeze(vii(1,4:5,2)));
% end
% 
% for ii = 2:Nv_f
%     im = squeeze(ims(:,:,:,ii));
%     imgs = rgb2gray(im);
%     
%     %     vii(ii,:,2) = vii(1,:,2);
%     if doplots
%         set(h1,'CData',im);
%         set(h2,'XData',squeeze(vii(ii,1:2,1)),'YData',squeeze(vii(ii,1:2,2)));
%         set(h3,'XData',squeeze(vii(ii,2:3,1)),'YData',squeeze(vii(ii,2:3,2)));
%         set(h4,'XData',squeeze(vii(ii,3:4,1)),'YData',squeeze(vii(ii,3:4,2)));
%         set(h5,'XData',squeeze(vii(ii,4:5,1)),'YData',squeeze(vii(ii,4:5,2)));
%         drawnow;
%     end
% end

% [~,vname,~] = fileparts(vpath);
% vo = VideoWriter(fullfile(userpath,'personal/Assignment',[vname '_ssd_out']),'Uncompressed AVI');
% vo.FrameRate = 10;  % Default 30
% open(vo);
% 
% for ii = 1:Nv_f
%     writeVideo(vo, double(squeeze(dnm(ii,:,:)./max(dnm(:)))));
% end
% close(vo);

    dnmx = dnm(:,33+[-8:8],33+[-8:8]);
    imagesc(squeeze(cat(3,dnmx(14,:,:),dnmx(15,:,:),dnmx(16,:,:),dnmx(17,:,:),dnmx(18,:,:))));aet;
    xlabel('Frame number');
%     ylabel('Vertical shift (pixels)');
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    title('SSD');
    colorbar;
    saveas(1,'SSD_GD_plot.png');

[~,vname,~] = fileparts(vpath);
vo = VideoWriter(fullfile(userpath,'personal/Assignment',[vname '_out']));
vo.FrameRate = v.FrameRate;  % Default 30
open(vo);

for ii = 1:Nv_f
    im = squeeze(ims(:,:,:,ii));
    vt = squeeze(p_rect(ii,:,:));
    M = [vt(1,1:2) vt(2,1:2);...
        vt(2,1:2) vt(3,1:2);...
        vt(3,1:2) vt(4,1:2);...
        vt(4,1:2) vt(1,1:2)];
    imo = insertShape(im,'Line',M,'LineWidth',5,'Color',[1 0 0]);
    writeVideo(vo, imo);
end
close(vo);

% end

