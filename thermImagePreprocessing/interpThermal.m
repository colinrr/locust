function interpThermal(matDir,interpDir,Tpath,geomf,Idx,dxdz,vmax,polyFile)
% interpThermal(matDir,interpDir,Tpath,geomf,Idx,dxdz,vmax,polyFile)
%       matDir = directory of raw .mat files from irbAsc2Mat
%       interpDir = directory to save new interpolated images
%       T         = path to table of image file names and indices -
%                   OPTIONALLY cell with a SECOND path to table containing
%                   plume masks. If any frames in Idx do not have
%                   associated masks, they will be generated via
%                   interpmask. This is useful when plumeTracker is run
%                   with some frames skipped (often necessary). This allows
%                   the two tables to be merged, and masks generated for 
%                   all frames. 
%                   (i.e. all frames, no masks + fewer frames, with masks)
%       geomf     = path to geometry file from mapPixels
%       Idx       = indices to run (OPTIONAL, defaults to all in T)
%       dxdz      = [OPTIONAL] [dx dz] Hard code the x,z spacing in interpolated image
%                   If scalar, makes dx=dz. By default, rounds to nearest
%                   0.1 of calculated value.
%       vmax      = [OPTIONAL] an INITIAL estimate for a maximum "realistic"
%                   object velocity, for use in smoothing the output
%                   masks. Any boundary pixels that "move" relative to the previous
%                   mask faster than this threshold are considered non-physical
%                   and removed. If left empty, will use max plume-top
%                   velocity * 2, if they are calculated. Otherwise, will
%                   leave empty and no smoothing occurs.
%       polyFile  = [OPTIONAL] file with manual polygons for clipping masks
%
% Note: It's important that the first plume masks be pretty accurate (ie
% actually capture the plume) for the mask/smoother interpolator to work
% properly.
%
% C Rowell, June 2018
fprintf('\n========= Interpolate Thermal =========\n')

narginchk(4,8)

if nargin<8
    polyFile = [];
end
if nargin<7
    vmax = [];
end
if nargin<6
    dxdz = [];
end
if isempty(dxdz)
    dxdz = 0;
else
    assert(or(isscalar(dxdz),numel(dxdz)==2),'Input "dxdz" must be scalar or 2-element vector.')
end

if nargin<5
    Idx = [];
end

if iscell(Tpath)
    load(Tpath{2}) % Load plumeTrackOutput (probably with subset of frames)
    Tmask = T.Mask;
    TmaskI = str2double(T.Properties.RowNames);
    Tpt = T; % PlumeTrack table
    
    if isempty(vmax)
        % Assuming plumeTrackCalcs is done, give a factor of 2 "safety" margin
        vmax = nanmax(T.v)*2;
    end

load(Tpath{1}) % Loads original table
elseif ischar(Tpath)
    load(Tpath); % Currently assumes masks are present if only one table entered
    Tmask = T.Mask;
    TmaskI = str2double(T.Properties.RowNames);
else
    error('Input path(s) to table(s) not recognized.')
end


% load(T)
load(geomf, 'geom') 
if ~isempty(polyFile)
    load(polyFile)
    use_poly = true;
else
    use_poly = false;
end

%% MERGE TABLES
interpVars = {'H','v','v_smooth','Ap'}; % Table vars to merge and interpolate from plumeTrack
killVars  = {'msec','X','Y'}; % Junk vars to remove from non-plumeTrack
% Also generate: new VidTime, interpFlag

% Get indices and cut table
if ~isempty(Idx)
    T = T(cellstr(string(Idx)),:);
end
Idx = double(string(T.Properties.RowNames));
N   = numel(Idx);

% Get rid of extra frame 
keepFrames = ismember(str2double(Tpt.Properties.RowNames)',Idx);
Tpt = Tpt(keepFrames,:); % Clear extra frames
interpFrames = ~ismember(Idx,str2double(Tpt.Properties.RowNames)');

% Clear junk vars
for jk = 1:length(killVars)
    if ismember(killVars{jk},T.Properties.VariableNames)
        T.(killVars{jk})=[];
    end
end

warning('off','MATLAB:interp1:NaNstrip')
for iv=1:length(interpVars)
    if ismember(interpVars{iv},Tpt.Properties.VariableNames)
        nanI = isnan(interp1(Tpt.Time,Tpt.(interpVars{iv}),T.Time));
        T.(interpVars{iv}) = interp1(Tpt.Time,Tpt.(interpVars{iv}),T.Time,'pchip');
        T.(interpVars{iv})(nanI) = NaN;
    end
end
% warning('on')

T.Time = T.Time-T.Time(1); % Reset time vector to zero
T.interpFlag = interpFrames; % Track which frames have interpolated values

%% START MAIN LOOP THROUGH FRAMES
count = 0;
fprintf('Source: \t%s\n',matDir)
fprintf('Destination:\t%s\n', interpDir)
fprintf('Number of images: %i\n\n',N)

F(N).Frame          = []; % Frame struct
Fp(N).mask          = []; % Meta-data struct
Fp(N).FileDateTime  = [];
Fp(N).ofile         = [];
Fp(N).xpg           = [];
Fp(N).zpg           = [];
Fp(N).idx           = [];
Fp(N).t             = [];
Fp(N).x             = [];
Fp(N).z             = [];
Fp(N).mm            = [];
Fp(N).mi            = [];

Fint = []; % Interpolation masks

% p=parpool(3);
dispstat('','init')
for kk = 1:N
    Fp(kk).idx = Idx(kk);
    Fp(kk).t   = T.Time(num2str(Fp(kk).idx));
    
    [Fp(kk).mm,Fp(kk).mi] =  ismember(Fp(kk).idx,TmaskI);
    if Fp(kk).mm
        Fp(kk).mask = logical(full(Tmask{Fp(kk).mi}));
    else
        Fp(kk).mask = [];
    end
%     load(fullfile(inputDir,T.File{idx}))
    Fp(kk).ofile = fullfile(interpDir,['int_' T.File{num2str(Fp(kk).idx)}]);
    
    if exist(Fp(kk).ofile,'file')
%         fprintf('%i: Re-gridded image already exists\n',Fp(kk).idx)
        dispstat(sprintf('%i/%i  %i: Interpolated image already exists\n',kk,N,Fp(kk).idx),'keepthis')
        
        load(Fp(kk).ofile, 'mask', 'dx', 'dz');
    else
        dispstat(sprintf('%i/%i  %i: Re-gridding image...\n',kk,N,Fp(kk).idx))
        F(kk) = load(fullfile(matDir,T.File{num2str(Fp(kk).idx)}),'Frame');
        [Fp(kk).xpg,Fp(kk).zpg]=meshgrid(1:size(F(kk).Frame,2),1:size(F(kk).Frame,1));
        [Fp(kk).x,Fp(kk).z]=px2m(Fp(kk).xpg,Fp(kk).zpg,geom);
    
        if and( use_poly, ~isempty(Fp(kk).mask) )
            % Apply manual polygons here
            Fp(kk).mask = applyPolygon(Fp(kk).idx,Fp(kk).mask,Polys);
        end

        [Frame,mask,gx,gz,dx,dz] = gridThermalFrame(F(kk).Frame,Fp(kk).mask,Fp(kk).x,Fp(kk).z,dxdz,[]);
        
        if and(~isempty(mask),exist('prevMask','var'))
            % Mask smoother
            if ~isempty(vmax)
                mask = smoothmask(mask,prevMask,prevt,Fp(kk).t,dx,dz,vmax);
            end

            if and(exist('prevMask','var'), ~isempty(Fint))
                dispstat(sprintf('   Writing %i interpolated mask(s): %s',...
                    length(Fint),sprintf('%i ',[Fint.idx])),'keepthis','keepthis')
                currMask = mask;
                intMasks = interpPTmask(prevMask,currMask,prevt,[Fint.t],Fp(kk).t);
                for ll=1:length(Fint)
                    % Do the mask interpolations
                    mask = intMasks(:,:,ll);
                    save(Fint(ll).ofile,'mask','-append')
                end
                Fint = [];
                mask = currMask; % yeap purrty important
            end
        elseif exist('prevMask','var') % Needs a previous mask to interpolate
            % Assemble list of masks to interpolate between this mask and the previous
            Fint = [Fint Fp(kk)]; 
        end
        
        % Interpolate missing mask
       
        if ~isempty(Fp(kk).ofile) && ~and(kk==N,isempty(mask))
            mask = sparse(mask);
            save(Fp(kk).ofile,'Frame','mask','gx','gz','dx','dz');
            mask = full(mask);
            T.File{num2str(Fp(kk).idx)} = ['int_' T.File{num2str(Fp(kk).idx)}];
            count = count+1;
        elseif and(kk==N,isempty(mask))
            warning('Last index requires interpolated mask - skipping this frame.')
            T(num2str(Fp(kk).idx),:) = []; % clear out of table
        end       
    end
    
    if ~isempty(mask)
        prevMask = mask;
        prevt    = Fp(kk).t;
        pdx = dx; % Just to check they are the same for now
        pdz = dz;
    end
    % Clear out vars
    F(kk).Frame          = [];
    Fp(kk).mask          = [];
    Fp(kk).FileDateTime  = [];
    Fp(kk).ofile         = [];
    Fp(kk).xpg           = [];
    Fp(kk).zpg           = [];
    Fp(kk).idx           = [];
    Fp(kk).x             = [];
    Fp(kk).z             = [];
    Fp(kk).mm            = [];
    Fp(kk).mi            = [];
    
end
% delete(p)
fprintf('Images re-gridded: %i/%i\n\n',count,N)
% Output frame Headers
interp_output_time = datestr(now,'dd-mm-yyyy HH:MM:SS');
save(fullfile(interpDir,'frameHeads.mat'),'T','interp_output_time')
end

function mask = applyPolygon(idx,mask,Polys)
% Do the thing!

for pp=1:length(Polys)
    if ismember(idx,Polys(pp).Idx)
        mask = mask.*Polys(pp).Mask;
    end
end
end

function mask = smoothmask(mask, prevMask, t1, t2, dx, dz, vmax)
% Smooth masks using a distance transform, timestamps, and velocity
% threshold

if issparse(prevMask)
    prevMask = full(prevMask);
end
if issparse(mask)
    mask = full(mask);
end
mask0=mask;

aspect = [dz dx 1];
D  = bwdistsc(edge(prevMask,'sobel'),aspect);
dMdt = D.*xor(prevMask,mask)./(t2-t1);

% dM = D.*mask;

% dMdt = dM(abs(dM)>0)/(t2-t1);

%% Method 1 - restrict expansion to any pixels below velocity threshold
%  '-> "safe" but weird masks can expand over time
% mask(dMdt>vmax) = 0;

%% Method 2 - eliminate any block containing pixels beyond velocity threshold
%  '-> "risky" because it may chop some true movement, but eliminates the
%  slow expansion problem
% dMmask = dMdt>0;
% cutmask = find(dMdt>vmax);
% CC = bwconncomp(dMmask);
% pixIdx = CC.PixelIdxList;
% 
% cutCC = cellfun(@(x) any(ismember(cutmask,x)),pixIdx);
% cutIdx = pixIdx(cutCC);
% cutIdx = unique(cat(1,cutIdx{:}));
% 
% mask(cutIdx) = prevMask(cutIdx);

%% Method 3 - Eliminate pixels growing faster than vmax
%  '-> made safer now by refining vmax estimate based on non-thresholded
%  blocks
dMmask = dMdt>0;
cutmask = find(dMdt>vmax);
CC = bwconncomp(dMmask);
pixIdx = CC.PixelIdxList;

cutCC = cellfun(@(x) any(ismember(cutmask,x)),pixIdx);
refIdx = pixIdx(~cutCC);
refIdx = unique(cat(1,refIdx{:}));

vmax2 = max(dMdt(refIdx)); % Update vmax using statistics of other moving blocks
if ~isempty(vmax2) && vmax2<vmax
    vmax = vmax2;
end

% Ensure 1 shape
if and( any(prevMask(:)), any(mask(:)) )
    mask(dMdt>vmax) = prevMask(dMdt>vmax);
    mask = biggestConnexComponent(mask);
end

end

function masks = interpPTmask(prevMask,mask,prevt,tint,t)
% use shape interpolation to generate masks between plumeTrack time steps
% prevMask, mask = previous and current plumeTrack masks
% prevt          = timestamp of previous mask
% tint           = timestamp(s) of missing masks to interpolate
% t              = timestamp of current mask

if ~islogical(prevMask)
    prevMask = logical(prevMask);
end
if ~islogical(mask)
    mask = logical(mask);
end

if issparse(prevMask)
    prevMask = full(mask);
end
if issparse(mask)
    mask = full(mask);
end

masks = interpmask([prevt t],cat(3,prevMask,mask),tint);

end
