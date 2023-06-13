function [D,opath] = getThermCube(inputDir,param,geomf,Idx,ROI,odir)
% [D,opath] = getThermCube(inputDir,param,geomf,Idx,ROI,odir)
% Uses SPATIALLY INTERPOLATED image frames to produce a Temperature data
% cube in x,z,t space by resampling in the time dimensions.
%
%       inputDir    = directory containing INTERPOLATED .mat frame files
%                       from plumeSpec1D - MUST ALL BE THE SAME SIZE!
%       param       = path to header file for inputDir frames
%       geom        = output geometry file from pixel mapping. Used mainly
%                       to reset vertical coordinate to vent-centered
%                       (assumes currently camera-centered).
%                       Leave empty to keep existing vertical coordinate 
%                       system, or enter a custom float to subtract from 
%                       vertical coordinates.
%       idx         = frame indices to gather [defaults to all]
%       ROI         = [1x4] array of pixel coordinates for thermal data gather
%                       (defaults to greatest mask extent)
%                       [x1 x2 z1 z2];
%       odir        = location to save thermal data cube. If empty, cube is
%                        not written, but is always ouput by the function.
%
%       D = output structure with fields:
    %       tname:  name of time vector from original header table. Usually 'Time'
    %       matDir: directory of individual frame files used to build dataCube
    %       refLLE: Lat,Lon,Elevation for frame reference pixel
    %       obsLLE: Lat,Lon,Elevation for camera
    %       targetLLE: Lat,Lon,Elevation for coordinate origin/frame plane
    %       z0ASL:  Reference height, m Above Sea Level, for z coordinate vector
    %       t_orig: time vector from original frames
    %       t0:     time vector starting value
    %       xp:     thermal profile pixel x vector
    %       zp:     thermal profile pixel z vector
    %       idx:    index vector from original header table
    %       imszOriginal:   original image sizes
    %       ROI:    pixel crop bounds from original images
    %       x:      thermal profile distance x array [m from center]
    %       z:      thermal profiles height z array [m above Z0]
    %       t:      regularly resampled time vector
    %       dx:     x grid spacing [m]
    %       dz:     z grid spacing [m]
    %       T:      time resampled cubic array Temperature values with i,j,k = z,x,t
    %       mask:   cubic logical array of mask values  i,j,k = z,x,t
    %
    %   These fields currently deprecated or left empty. Can be uncommented
    %   below:
    %       dT:     Gradient of T along t (dT/dt)
    %       I:      Converted to effective at-sensor radiance (simple Stefan
    %               Bolztmann, I = epsilon*sigma*T.^4, with epsilon=1)
    %       dI:     Gradient of I along t (dI/dt)
    %       Tmax:   Horizontal (across x) maximum of T
    %       Tint:   Horizontally (across x) integrated T
    %       Iint:   Horizontally (across x) integrated I
    %       dAint:  Horizontally (across x) integrated dT
    %       dIint:  Horizontally (across x) integrated dI
%       opath = full path to save file
%

%       
% C Rowell 2018
%
fprintf('\n========= getThermCube =========\n')
clear textprogressbar
%% Parse input
if nargin<6
    odir = [];
end
if nargin<5
    ROI = [];
end
if nargin<4
    Idx = [];
end    
if nargin<3
    geomf = [];
end

defaultSmoothLengthSeconds = 4.5; % Default smoothing length frame MASKS (# of frames)
%% Load parameter files and set defaults
disp('Loading frame parameters and geometry...')

if isempty(geomf)
    Z0 = 0;
elseif ischar(geomf)
    load(geomf,'geom')
    Z0 = geom.Ztarg-geom.Z0; % Assumes camera-centered z coordinate in input frames
elseif and(isnumeric(geomf),numel(geomf)==1)
    Z0 = geomf;
else
    error('Geom input not recognized.')
end
% error('flargh: Dev. note: Make sure geom information and Z0 are properly implemented next time you run this')
% Probably also requires properly fixing mapPixels and geom output

load(param)

D.tname = 'Time';

if ~isempty(Idx)
    T = T(cellstr(string(Idx)),:);
end
Idx = T.Properties.RowNames;

if isempty(ROI)
    ROI=cellfun(@(x) [find(sum(x,1),1,'first') find(sum(x,1),1,'last') find(sum(x,2),1,'first') find(sum(x,2),1,'last')], T.Mask,'UniformOutput',false);
    ROI =cell2mat(ROI(~cellfun(@isempty,ROI)));
    ROI = [min(ROI(:,1)) max(ROI(:,2)) min(ROI(:,3)) max(ROI(:,4))];
end

%%
sig = 5.67e-8;

D.matDir = inputDir;

if exist('geom','var')
%     D.refLLE = geom.refLLE;
%     D.obsLLE = geom.obsLLE;
%     D.targetLLE = geom.targetLLE;
    D.geom = geom;
end
D.z0ASL = geom.Ztarg;
% Get reference name (string idx for later use)
D.t_orig = T.(D.tname)(Idx);
D.t0_local = T.Timestamp(Idx{1});
D.t0 = min(D.t_orig);



D.xp = ROI(1):ROI(2);
D.zp = ROI(3):ROI(4);

Mz       = numel(D.zp);
Nx       = numel(D.xp);
Ot       = numel(Idx); %size(T,1);

A     = zeros([Mz Nx Ot]);
Amask = zeros([Mz Nx Ot]);
D.idx   = Idx;

%% Main loop to collect image data
fprintf('Retrieving temperatures for %i frames...\n',Ot)
textprogressbar('   ... ')
for kk = 1:Ot
    idx = Idx{kk};
    
    load(fullfile(inputDir,T.File{idx}))
    if kk==1
        D.imszOriginal = size(Frame);
        D.ROI = ROI;
        

        
    end
    
    % ------ Pull frames ------
    A(:,:,kk) = Frame(D.zp,D.xp);

    if ~isempty(mask)
        Amask(:,:,kk) = mask(D.zp,D.xp);
    else
        Amask(:,:,kk) = zeros(numel(D.zp),numel(D.xp));
    end
    
    textprogressbar(kk/Ot*100)
end
textprogressbar(' -> Done')

%% Spatial vectors and coordinate flip
D.t = D.t_orig - min(D.t_orig); % Reset time to zero
D.x = gx(D.xp)'; % Using interpolated vectors
D.z = gz(D.zp)';
D.z = D.z - Z0;
D.dx = dx;
D.dz = dz;

% Change from image pixel coordinates to plume coordinates
D.z = flipud(D.z);
D.dz = -D.dz;
A = flipud(A);
Amask = logical(flipud(Amask));


%% Reshape, resample, smooth
disp('Reshaping and resampling data cube...')

% Data array
A = permute(A,[3 1 2]);
A = reshape(A,[Ot Mz*Nx]);
[A,D.t] = resample(A,D.t,'pchip');
A = reshape(A,[Ot Mz Nx]);
D.T = single(permute(A,[2 3 1])); % Probably don't need double precision...?
clear A

% Mask
Amask = permute(Amask,[3 1 2]);
Amask = reshape(Amask,[Ot Mz*Nx]);
Amask = interp1((D.t_orig-min(D.t_orig)),single(Amask),D.t,'nearest',0);
Amask = reshape(Amask,[Ot Mz Nx]);
D.mask = logical(permute(Amask,[2 3 1]));
clear Amask

%% Mask smoothing
disp('Smoothing masks...')
area = squeeze(sum(D.mask(:,:,1:end-1),[1 2]));
sLengths = [3:2:39]; % Range of smoothing lengths to check

% Get mask area (# pixels) for each frame as a test proxy for mask smoothness
As = zeros(length(area),length(sLengths));
SAD = zeros(length(sLengths),1);
for ii=1:length(sLengths)
    As(:,ii)  = round(smooth(area,sLengths(ii))); % Approximation to smoothed mask areas
    SAD(ii) = sum(abs(As(:,ii) - area)); % Summed absolute differences (smoothed area - raw area)
end
dSAD = 1-(SAD(1:end-1)-min(SAD))./(SAD(2:end)-min(SAD));

% Get optimal smoothing length
Isad = find(dSAD<0.1,1,'first');
if isempty(Isad)
    D.maskSmoothLength = 2*round(defaultSmoothLengthSeconds/mean(diff(D.t))/2)+1; % Nearest odd number
    [~,Isad] = closest(D.maskSmoothLength,sLengths);
else
    if Isad==1
        Isad = Isad+1; % Smoothing length of 3 should always be rejected
    end
    D.maskSmoothLength = sLengths(Isad);
end

% Temporary plot to make sure this part goes "smoothly"...
figure
subplot(2,1,1)
plot(sLengths(1:end-1),dSAD)
hold on
plot(D.maskSmoothLength,dSAD(Isad),'or')
xlabel('Smoothing length')
ylabel('delta(Summed Absolute Differences)')
title('Mask smoothing results')
subplot(2,1,2)
plot([area,As(:,Isad)])
xlabel('Frame')
ylabel('Mask area (# pixels)')
legend('Raw masks',sprintf('Smoothed, N = %i, dSAD = %.2f',D.maskSmoothLength,dSAD(Isad)))

% Run smoothing
D.mask = smooth3(double(D.mask),'box',[3 3 D.maskSmoothLength])>0.5;
% Conncomp to eliminate any small separated chunks
for ii=1:size(D.mask,3)
    if sum(D.mask(:,:,ii),[1 2])>0
        D.mask(:,:,ii) = biggestConnexComponent(D.mask(:,:,ii));
    end
end
% POSSIBLY use regionprops to fill any holes

%% Calculate gradients, maxima, normalizations, etc
% disp('Calculating secondary data fields...')

% Un-comment items below to add them into the routine if desired



% D.Amask = sparse(D.Amask);
% D.Tmax = squeeze(max(D.T,[],2)); % Max temp across x
% D.Tmax = [];

% Fast integration approach, requires interpolated images-----
% Raw integration
% D.Tint = squeeze(trapz(D.x,D.T,2)); % Temp integrated across x
% D.Tint = [];

% ----- Temperature Gradient (+ x integration) ------
% [~,~,D.dT] = gradient(D.T,D.x,D.z,D.t_orig);
% D.dTint = squeeze(trapz(D.x,D.dT,2));

% ---- Using simple diff -----
% dA = diff(D.A,1,3);
% dAint = squeeze(trapz(D.x,dA,2));
% -------------------------------

% ---- Try converting to effective at-sensor radiance ------
% D.I = sig*D.T.^4;
% D.I = [];
% [~,~,D.dI] = gradient(D.I,D.x,D.z,D.t_orig); % Gradient of radiance

% Horizontal integrations
% D.Iint = squeeze(trapz(D.x,D.I,2)); % Integrated radiance
% D.Iint = [];
% D.dIint = squeeze(trapz(D.x,D.dI,2)); % Integrated gradient of radiance
% D.Aint = D.Aint-D.Aint(:,1); % Subtract reference integration


if and(~isempty(odir),exist(odir,'dir'))
    ofile = sprintf('thermStats_%s_z%i_x%i_t%i',datestr(now,'YYYY-mm-dd'),Mz,Nx,Ot);
    opath = fullfile(odir,ofile);
    fprintf('Writing:  %s\n',opath)
    save(opath,'D','-v7.3')
    opath = [opath '.mat'];
else
    opath = [];
end
disp('.....DONE!')
