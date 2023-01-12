% Plot a basic illustration of the build-up of the reconstruction in
% (unfiltered) back-projection.
%
% Samuli Siltanen Sep 2020

%% Preliminaries

% Parameters for controlling the appearance of the plot
lwidth = .5;
detwidth = .8;
raywidth = .5;
datalinewidth = 1.5;
RecoBoundaryWidth = .7;
RecoBoundaryColor = [23 91 184]/255;
GuideLineColor    = .4*[1 1 1];
projamp = 400; % Amplitude of projection data in the plot
projoff = 15;  % Offset of projection data in the plot
scmN   = 30;
gammacorr = .8; % Gamma correction for showing the Shepp-Logan phantom
purple   = [186 0 255]/255;
yellow   = [246 255 0]/255;
detcolor = [70 173 211]/255;

% Build phantom
Nph = 1024;
ph  = double(imread('data/banana.jpg','jpg'));
ph  = max(0,ph);
ph  = ph/max(ph(:));

% Select rotation angles for the image
Nang = 400;
rotangs   = [0:(Nang-1)]/(Nang-1)*pi;  % Angles in radians
rotangs_d = [0:(Nang-1)]/(Nang-1)*180; % Angles in degrees
% Compute and save Radon transform
[sinogram,t] = radon(ph,-90+rotangs_d);
MAXsino      = max(max(sinogram));

% Define a couple of guidelines indicating the direction of X-rays by
% specifying starting points and end points
Nrays = 2;
startX = -.5*Nph*ones(1,Nrays);
endX   = 1.3*Nph*ones(1,Nrays);
startY = 3*Nph/2 + t(1) + (t(end)-t(1))*[0:(Nrays-1)]/(Nrays-1);
endY   = startY;

% Open video file
videofilename = 'banana_BP';
videotype = 'MPEG-4';
v1 = VideoWriter(videofilename,videotype);
v1.Quality = 95;
open(v1);

%% Create frames

% Loop over rotation angles
for jjj = 1:Nang
    
    % Create plot window
    figure(1)
    clf
    
    % Draw back-projection using the tomographic data collected up to the
    % current angle
    tmpsino = zeros(size(sinogram));
    tmpsino(:,1:jjj) = sinogram(:,1:jjj);
    BP = iradon(tmpsino,rotangs_d-90,'linear','none');
    BP = BP(2:(end-1),2:(end-1));
    BP = max(0,BP);
    BP = BP/max(BP(:));
    BP_plot = [zeros(Nph);BP;zeros(Nph)];
    BP_plot = [zeros(3*Nph,Nph),BP_plot,zeros(3*Nph,Nph)];
    imagesc(BP_plot.^gammacorr,[0,1])
    colormap gray
    hold on
    
    % Construct rotation matrix
    th = rotangs(jjj);
    rmat = [[cos(-th),-sin(-th)];[sin(-th),cos(-th)]];
    
    % Rotate X-rays
    tmp1 = [startX-Nph/2;startY-3*Nph/2];
    tmp2 = [endX-Nph/2;endY-3*Nph/2];
    tmp1 = rmat*tmp1;
    tmp2 = rmat*tmp2;
    startXr = tmp1(1,:)+3*Nph/2;
    startYr = tmp1(2,:)+3*Nph/2;
    endXr = tmp2(1,:)+3*Nph/2;
    endYr = tmp2(2,:)+3*Nph/2;
    
    % Plot the guidelines
    for iii = 1:Nrays
        plo1 = plot([startXr(iii), endXr(iii)],[startYr(iii), endYr(iii)],'k','linewidth',raywidth);
        set(plo1,'color',GuideLineColor);
    end
    
    % Plot detector (using guideline color to de-emphasize it)
    plo3 = plot([endXr(1) endXr(end)],[endYr(1) endYr(end)],'k','linewidth',detwidth);
    set(plo3,'color',GuideLineColor);
    
    % Plot projection data on top of detector
    projdata = sinogram(:,jjj);
    paravec  = [endXr(end);endYr(end)]-[endXr(1);endYr(1)];
    perpunitvec = [paravec(2);-paravec(1)];
    perpunitvec = perpunitvec/norm(perpunitvec);
    tmp = linspace(0,1,length(projdata));
    plotpoints = zeros(2,length(projdata));
    for iii = 1:length(projdata)
        plotpoints(:,iii) = [endXr(1);endYr(1)]+tmp(iii)*paravec+(projoff+projamp*(projdata(iii)/MAXsino))*perpunitvec;
    end
    plo4 = plot(plotpoints(1,:),plotpoints(2,:),'linewidth',datalinewidth);
    set(plo4,'color',yellow);
    
    % Plot boundaries of the imaging domain
    p7 = plot([Nph 2*Nph 2*Nph Nph Nph],[Nph,Nph 2*Nph 2*Nph Nph],'linewidth',RecoBoundaryWidth);
    set(p7,'color',RecoBoundaryColor);
    
    % Axis settings
    axis([-Nph 3*Nph -Nph 3*Nph])
    axis equal
    axis off
    
    % Initial image 
    im1 = print('-r600','-RGBImage');
    [row1,col1] = size(im1(:,:,1));
    
    % Indices of sides of the black square
    bs1_left   = min(find(im1(round(row1/2),:,1)==0));
    bs1_right  = max(find(im1(round(row1/2),:,1)==0));
    bs1_top    = min(find(im1(:,round(col1/2),1)==0));
    bs1_bottom = max(find(im1(:,round(col1/2),1)==0));
    
    % Crop the image
    startrow = round(bs1_top+.15*Nph);
    endrow   = round(bs1_bottom-.65*Nph);
    startcol = round(bs1_left+.2*Nph);
    endcol   = round(bs1_right-.2*Nph);
    im1 = im1(startrow:endrow,startcol:endcol,:);
        
    % Adjust image height to 1080
    im2 = uint8(zeros(1080,1920));
    im3 = imresize(im1, [1080 NaN]);
    [row3, col3, tmp3] = size(im3);
    im2(:,round((1920-col3)/2)+[1:col3],1) = im3(:,:,1);
    im2(:,round((1920-col3)/2)+[1:col3],2) = im3(:,:,2);
    im2(:,round((1920-col3)/2)+[1:col3],3) = im3(:,:,3);
    
    % Add frame to video
    writeVideo(v1,im2);
    
    % Monitor the run
    if mod(jjj,10)==0
        disp([jjj,Nang])
    end
    
end



%% Add freeze frame

% Duplicate the last frame
for jjj = (Nang+1):4*Nang
    
    % Add frame to video
    writeVideo(v1,im2);
    
    % Monitor the run
    if mod(jjj,10)==0
        disp([jjj,4*Nang])
    end
end

close(v1);
