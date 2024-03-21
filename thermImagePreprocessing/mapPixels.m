function [geom,savefile] = mapPixels(obs,target,ref,ij,hfov,vfov,imsz,regParams,savepath)
% [geom,savefile] = mapPixels(obs,target,ref,ij,hfov,vfov,imsz,regParams,savepath)
% Function to create a mapping for pixel position to distance. Assumes
% image to be in a 2D plane.
% IN:   obs  = [lat, lon, elev] coords of the camera observing position
%       target = [lat, lon, elev] point through which to fix the projected
%               image plane
%       ref  = [lat, lon, elev] coords of a reference point in the image
%       ij   = [i j] pixel indices of same reference point in image
%       hfov = horizontal field of view (degrees)
%       vfov = vertical field of view (degrees)
%       imsz = num. pixels [y,x]
%       regParams = [optional] Struct or path to struct to use results of image 
%                  registration for pixel cropping/image bounds.
%       savepath = [optional] directory to save output geometry structure to .mat.
% 
% OUT:  
%       geom = a struct containing geometric information - distances and
%           angles from the observation point, frame geometry, etc.
%
%       savefile = full path with filename for saved file. Only valid if
%            savepath is entered as input.

% Needs to use a reference point in the image and camera location to get 
% inclination angle of camera, then calculates a spatial mapping based on FOV
%
% C Rowell, July 2018
% Based on Harris, 2013, Thermal Remote Sensing of Active Volcanoes: 
% A User's Manual, Chapter 9

fprintf('\n========= Map Pixels to Local Coordinates =========\n')

if nargin<9
    savepath = [];
end
if nargin<8
    regParams = [];
end
if exist('regParams','var')
    if ~isempty(regParams)
        use_registration = true;
    else
        use_registration = false;
    end
else
    use_registration = false;
end

%%
% Fix vertical plane of image above vent for now. In general, this is a
% reference point in the desired plane of interest, while "ref" is simply
% a known fixed point contained in the image at arbitrary distance.
% vent = [-15.786744, -71.855919, 5911];

spheroid = referenceEllipsoid('WGS 84');

% Thetas can be elevation angles, phis can be azimuth angles
[az,thetaR_prime,Dr] = geodetic2aer( ...
    ref(1),ref(2),ref(3),obs(1),obs(2),obs(3),spheroid);
[azV, thetaV_prime, Dv] = geodetic2aer( ...
    target(1),target(2),target(3),obs(1),obs(2),obs(3),spheroid);

% The above elevation angles given "apparent" elevation angle from camera
% point. Corrections below are needed for "true" elevation angles from
% camera axis.

% Get total vertical and horizontal fields of view, find distance to view
% center that corresponds to reference point. May later want to project
% this back/forward to estimated vent location or something similar.
% Recover camera elevation angle from landmark

ny   = imsz(1);
nx   = imsz(2);
vifov=vfov/ny; % Individual pixel fields of view
hifov=hfov/nx;

% Needs to account for odd/even image size - if not here then for sure in px2m
ny0 = ny/2;
nx0 = nx/2;

% !!! Need to carefully consider pixel centers vs edges here !!!
dTheta_ref      = (ij(1)-(ny0+0.5))*vifov;                  % true vertical angular distance between reference point and camera axis
dPhi_ref        = (ij(2)-(nx0+0.5))*hifov;                  % true angular distance between ref point and camera center along camera axis
thetaR          = atand(tand(thetaR_prime)/cosd(dPhi_ref)); % TRUE elevation angle of ref corrected for camera orientation
Theta0          = thetaR + dTheta_ref;                      % true Camera inclination angle (to image centerline: edge b/w pixels)
dPhi_ref_prime  = atand(cosd(Theta0)*tand(dPhi_ref));       % Azimuthal angular distance, accounting for camera tilt

Phi0            = az - dPhi_ref_prime;                      % Azimuth to camera centerline
dPhi_targ_prime = azV - Phi0;                               % Azimuthal angle between target point and viewpoint centerline.
dPhi_targ       = atand(tand(dPhi_targ_prime)/cosd(Theta0));% TRUE horizontal angular distance between target and centerline in camera reference frame
thetaV          = atand(tand(thetaV_prime)/cosd(dPhi_targ));   % TRUE elevation angle to target

Dnorm      = Dv*cosd(dPhi_targ)*cosd(thetaV);               % Horizontal distance to camera center pixel,equal to xRange unless target plane is tilted from vertical
Dlos       = Dnorm/cosd(Theta0);                            % Distance to image center point in target plane

VFOV       = Dnorm*(tand(Theta0+vfov/2)-tand(Theta0-vfov/2));% Total FOV extent in target plane, meters
Lp          = Dnorm*(tand(Theta0+vifov/2)-tand(Theta0-vifov/2)); % Approx center pixel dimension (if it were centered in the image)


% DEPRECATED THINGS---------
% dPhi_rv    = azV - az; % Angular distance between reference point and target point.
% ThetaNorm  = Theta0 - targetPlaneTilt; % Camera tilt relative to target plane;

% Dlos       = Dv*cosd(dPhi_targ)*cosd(thetaT)/cosd(Theta0);  % Distance to image center point in target plane
% Dlos       = Dv*cosd(thetaV)/cosd(Theta0);

% xRange     = Dlos*cosd(Theta0); 
% xRange     = Dnorm/cosd(targetPlaneTilt); % Min horizontal distance to target plane intersection

% y0         = Dlos*sind(Theta0); % Vertical distance above camera
% HFOV0      = [];
% HFOV1      = [];

% Size range of pixels
% Yn = xRange*( tand(Theta0 + (ny0 - (0:ny-1)')*vifov) - tand(Theta0 + (ny0 - (1:ny)')*vifov) ) ;
% Yth = Theta0 + (ny0-(1:ny))*vifov;
           
% Test the vis
% [icell,jcell] = meshgrid(1:nx,1:ny);
% [inode,jnode] = meshgrid(0.5:nx+0.5,0.5:ny+0.5);
% Xc = px2x(icell,jcell); Yc = px2y(jcell);
% Xn = px2x(inode,jnode); Yn = px2y(jnode);
% [x,y] = xzmap(inode,jnode);
% ---------------------------

geom.obsLLE         = obs;
geom.targetLLE      = target;
geom.refLLE         = ref;
geom.ref_pix_ij     = ij;
geom.target_pix_ij  = round([-(thetaV-Theta0)/vifov+ny0+0.5 dPhi_targ/hifov+nx0+0.5]);
geom.nx0            = nx0;
geom.ny0            = ny0;
geom.im_size        = imsz;
geom.vfov_deg       = vfov;
geom.hfov_deg       = hfov;
geom.vifov_px_deg   = vifov;
geom.hifov_px_deg   = hifov;
geom.VFOV_m         = VFOV;

geom.center_pixSz_m = Lp;
geom.center_azim    = Phi0;
geom.D_los          = Dlos;
geom.X_distance     = Dnorm;
geom.Elev_angle     = Theta0;
geom.Ztarg          = target(3); % Camera-centered coordinates Z origin (m ASL)
geom.Z0             = obs(3); % Target-centered coordinates Z origin (m ASL)

% geom.X_distance   = xRange;
% geom.h_distance   = xRange;
% geom.HFOV0_m        = HFOV0;
% geom.HFOV1_m        = HFOV1;% geom.Tilt_angle=ThetaNorm; % Usage needs updating

% geom.cam_LLE      = obs;

% Apply registration correction??
if use_registration
    disp('Applying geometry correction for image registration.')
    geom = registerPixelMap(geom,regParams);
end

disp('Camera position and orientation:')
fprintf('Slant distance =\t%f m\n',geom.D_los)
fprintf('Hor. distance  =\t%f m\n',geom.X_distance)
fprintf('Elev. angle    =\t%f degrees\n',geom.Elev_angle) 
fprintf('VFOV Extent    =\t%f m\n',geom.VFOV_m)
fprintf('Lp             =\t%f m\n',geom.center_pixSz_m)
% Image center is W/2+0.5, H/2+0.5 as plotted in matlab (pixel centered
% coords)

%% Test figure to verify the math

plotImageMap(geom)

% Save?
if ~isempty(savepath)
    if use_registration
        savefile = fullfile(savepath,sprintf('geometry_%s_registered.mat',datestr(now,'YYYY-mm-dd')));
    else
        savefile = fullfile(savepath,sprintf('geometry_%s.mat',datestr(now,'YYYY-mm-dd')));
    end
    fprintf('Saving geometry file:\n\t%s\n',savefile)
    save(savefile,'geom')
else
    savefile = [];
end

end