function lp=image2dem(frame,mask,geom,demfile,roi,matdir,frames)
%
%
%
%
%

% load(fullfile(outputDir,'output_params.mat'));
[A,X,Y] = loadDEM(demfile,roi);

% Obs pt
obs = geom.cam_LLE;
[obsN,obsE,~] = ell2utmWGS84(obs(1), obs(2));

% Get image frame coords
fj = [0.5 geom.im_size(2)+0.5; 0.5 geom.im_size(2)];
fi = [0.5 0.5; geom.im_size(1) geom.im_size(1)];
[fx,fz] = px2m(fj,fi);
[fE,fN,fZ] = xz2utm(fx,fz,geom);


% [flat,flon,fz] = px2geo(fx,fy);
% [fN,fE,~]     = ell2utmWGS84(flat,flon);
% [obsN,obsE,~] = ell2utmWGS84(obsLLE(1), obsLLE(2));
% obsZ = interp2(X,Y,A,obsE,obsN);

% Get wireframe for image pixels
dpx = 1;
px = 1:dpx:geom.im_size(2);
pz = 1:dpx:geom.im_size(1);
[px,pz] = meshgrid(px,pz);

% [plat,plon,pz] = px2geo(px,py);
% [pN,pE,~]     = ell2utmWGS84(plat,plon);
[x,z] = px2m(px,pz);
[E,N,Z] = xz2utm(x,z,geom);

% Plume outline from geodetic
% frame_num = 300;
% outline = dat.Outline{frame_num};
% [olat,olon,oz] = px2geo(outline(:,2),outline(:,1));
% [oN,oE,~]     = ell2utmWGS84(olat,olon);

% Plume outline from x,z
% [pE,pN,pZ] = xz2utm(x,z,geom);




if numel(frames)>1
    vid_mode = true;
    disp('Video output mode.')
    vidObj = VideoWriter(fullfile(matdir, 'PTresults/demPlume3'),'Motion JPEG AVI');
    vidObj.FrameRate = 20;
    vidObj.Quality = 100;
    open(vidObj);
else
    vid_mode = false;
end
E0 = min(X(:));
N0 = min(Y(:));
fig=figure;
set(fig, 'Position', [100 100 1*geom.im_size(2) 1*geom.im_size(1)])
set(gcf, 'PaperPositionMode', 'auto');

% view_param = [1 0.5 0.6];
% view_param = [1 0.8 0.25];
view_param = [0.2,1,0.2];

for fn = frames
    Frame=loadImg(fullfile(matdir,dat.File{fn})); % Go by frame number
    mask = logical(dat.Mask{fn});
    Frame(~mask)=NaN;
%% Plot up some jazz

    surf(X-E0,Y-N0,A,'EdgeAlpha',0,'FaceAlpha',1)
    daspect([1 1 1])
    camlight('left')
    material dull
    colormap(gray)
    axis tight
    view(view_param)
%     zl=zlim; zlim([zl(1) max(fZ(:))])
    zl=zlim; zlim([zl(1) 6500])
    cax = caxis; caxis(cax);
    xlabel('Rel Easting (m)')
    ylabel('Rel Northing (m)')
    hold on
%     plot3(obsE-E0,obsN-N0,obs(3),'y^','LineWidth',2)

    a=surf(fE-E0,fN-N0,fZ,'FaceAlpha',0.2);
    % freezeColors

    % Plot view lines
%     lx = [obsE*[1 1 1 1];fE(:)'];
%     ly = [obsN*[1 1 1 1];fN(:)'];
%     lz = [obs(3)*[1 1 1 1];fZ(:)'];
%     plot3(lx-E0,ly-N0,lz,'Color',[0.3 0.3 1])
    
    % figure
    axa = gca;
    axp = get(gca,'Position');
    axlims = axis;
    axb = axes('Position',axp);
    b = surf(E-E0,N-N0,Z,Frame(1:dpx:end,1:dpx:end),'EdgeAlpha',0);
    axb = gca;
    daspect([1 1 1])
    set(gca,'Color','none','XTickLabel',[],'YTickLabel',[],'YTickLabel',[])
    view(view_param)
    axis tight
    zlim([zl(1) max(fZ(:))])
    axis(axlims)
    grid off
    lp = linkprop([axa axb],'view');
    colormap(axb,hot)
    caxis([240 360])
    
    if vid_mode
        F = getframe(fig);
        writeVideo(vidObj,F);
        clf
    end
end
% a=surf(pE-E0,pN-N0,pZ,'FaceAlpha',0.2,'EdgeAlpha',0.7);
% b=plot3(oE-E0,oN-N0,oz,'r');
% a=surf(pE-E0,pN-N0,pz,Frame(1:dpx:end,1:dpx:end),'EdgeAlpha',0);
% a=scatter3(pE(:)-E0,pN(:)-N0,pz(:),10,Fr(:));
end

function [E,N,Z] = xz2utm(x,z,geom)
%     [obsN,obsE,~] = ell2utmWGS84(geom.cam_LLE(1), geom.cam_LLE(2));
%     obsZ = geom.cam_LLE(3);
    az = geom.center_azim-180; % Plane azimuth
%     if az<0; az = 360-az; end
    % Reference transformation equals coords of image center, elevation of
    % camera
    spheroid = referenceEllipsoid('WGS 84');
    [rlat,rlon,H] = aer2geodetic(geom.center_azim,0,geom.X_distance,...
                geom.cam_LLE(1),geom.cam_LLE(2),geom.cam_LLE(3),spheroid);
    [rN,rE,~]     = ell2utmWGS84(rlat,rlon);
    
    dTheta = 180 - az; % Rot. angle - assume starting with image plane azimuth at 180 (x,z plane)
    
    % Build transformation matrix
    T = diag([cosd(dTheta) cosd(dTheta) 1 1]);
    T(1,2) = -sind(dTheta);
    T(2,1) = sind(dTheta);
    T(1:3,4) = [rE,rN,H]';
    
    szx = size(x);
    
    C = [x(:)'; zeros([1 length(x(:))]); z(:)'; ones([1 length(x(:))])];
    
    Cp = T*C;
    E = reshape(Cp(1,:)',szx);
    N = reshape(Cp(2,:)',szx);
    Z = reshape(Cp(3,:)',szx);
end
