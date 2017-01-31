function vopath = videoTrack(vpath,coords)
% videoTrack Track rectangle in video file
%
% Given a video with moderate shakiness and coordinates of a rectangle in
% first frame, output a video with the rectangular section tracked over the
% course of the video.
%
% INPUT:
%
%   vpath:  "" full path to video file
%
%   coords: [4] rectangle coordinates specified by [x y w h] where x and y
%   are horizontal and vertical pixel offsets from the top left corner of
%   the video and w and h are the width and height of the rectangle, all
%   for frame 1 of the video file provided
%
% OUTPUT:
%
%   vopath: "" full path to output video file
%
% SEE:
%
% charles dunn
% 

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

No = 64;%buffer region to avoid objects going out of frame
Nn = 64;%number of offsets in each direction to check
Ns = 2;%size of local minimum check in each direction
Nt = .1;%angle resolution in degrees for rotation search

xst = 0;
yst = 0;

xs = zeros(1,Nv_f);
ys = zeros(1,Nv_f);
theta = zeros(1,Nv_f);
%iterate over frames
for ff = 1:Nv_f
    %calculate metrics
    [gradm, grada] = imgradient(squeeze(imgsf(No:end-No,No:end-No,ff)));
    grada = grada(gradm(:)>64);
    grada = mod(grada+90,180)-90;
    gradh = hist(grada,-90:Nt:90);
    gradh = gradh./sum(gradh);
    
    if(ff>1)
        %start search for SSD minimum at previous shift
        ssd = nan(65,65);
        ssd(yst+Nn/2+1,xst+Nn/2+1) = sum(colvec((...
            gradm_prev(1-min(0,yst):size(gradm_prev,1)+min(0,-yst),1-min(0,xst):size(gradm_prev,2)+min(0,-xst)) -...
            gradm(1-min(0,-yst):size(gradm,1)+min(0,yst),1-min(0,-xst):size(gradm,2)+min(0,xst))...
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
                        continue;%avoid repeated computation
                    end
                    if(isnan(ssd(ii+Nn/2+1,jj+Nn/2+1)))
                        ssd(ii+Nn/2+1,jj+Nn/2+1) = sum(colvec((...
                            gradm_prev(1-min(0,ii):size(gradm_prev,1)+min(0,-ii),1-min(0,jj):size(gradm_prev,2)+min(0,-jj)) -...
                            gradm(1-min(0,-ii):size(gradm,1)+min(0,ii),1-min(0,-jj):size(gradm,2)+min(0,jj))...
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
%         dnm(ff,:,:) = ssd;
        
%         imagesc(-8:8,-8:8,ssd(33+[-8:8],33+[-8:8]));aet;
%         xlabel('Horizontal shift (pixels)');
%         ylabel('Vertical shift (pixels)');
%         title('SSD');
%         colorbar;
%         saveas(1,'SSD_GD_plot2.png');
        
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
        
        %determine rotation between frames
        theta_shift = (-5/Nt:5/Nt);
        gradh_dot = zeros(1,length(theta_shift));
        for ii = 1:length(theta_shift)
            gradh_dot(ii) = sum(gradh_prev.*circshift(gradh,theta_shift(ii)));
        end
        [~,mind] = max(gradh_dot);
        
        theta(ff) = theta_shift(mind)*Nt;
        
        %apply rotation to rectangle coordinates
        R = [cosd(theta(ff)) -sind(theta(ff));sind(theta(ff)) cosd(theta(ff))];
        for ii = 1:4
            p_rect(ff,ii,1:2) = R*(squeeze(p_rect(ff,ii,1:2)) - [Nv_h/2; Nv_w/2]) + [Nv_h/2; Nv_w/2];
        end
        
    end
    if(ff~=Nv_f)
        gradm_prev = gradm;
        gradh_prev = gradh;
    end
end

% [~,vname,~] = fileparts(vpath);
% vo = VideoWriter(fullfile(userpath,'personal/Assignment',[vname '_ssd_out']),'Uncompressed AVI');
% vo.FrameRate = 10;  % Default 30
% open(vo);
%
% for ii = 1:Nv_f
%     writeVideo(vo, double(squeeze(dnm(ii,:,:)./max(dnm(:)))));
% end
% close(vo);

%write out video with tracked rectangle superimposed
[~,vname,~] = fileparts(vpath);
vopath = fullfile(userpath,'personal/Assignment',[vname '_out']);
vo = VideoWriter(vopath);
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

end

