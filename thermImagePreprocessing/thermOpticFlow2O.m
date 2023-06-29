function [V,vfile] = thermOpticFlow2O(thermCube,opticParams,ofile)
% [V] = thermOpticFlow(thermCube,opticParams,ofile)
% Uses Sun et al., 2010 optical flow package to calculate 2D velocity
% fields in thermal imagery. For end points (1st and last images), uses a
% forward/backward difference, respectively (ie uv(1) = f(Frame(1), Frame(2)).
% For all middle frames, uses a central difference (ie uv(i) = f(Frame(i-1), Frame(i+1)).
% This scheme gives a slightly more stable result.
% 
% INPUT:  thermCube = path to data cube struct from getThermCube
%         opticParams   = struct o fparamaters for any pre-processing of 
%                     images priorvto running optic flow
%             FIELDS:
%               Sub: integer vector of frame (time) subscripts in thermCube. 
%                   (NOT D.idx). Leave empty [] to take all in data cube
%               method: opticFlow method for Sun algorithm. Default
%                   is 'classic+nl-fastp'. See 'ijcv_flow_code/estimate_flow_demo'
%               maskPad: for efficient computation, frames are cropped to
%                   mask boundaries, plus a pad measured in pixels on all
%                   sides.
%               method = see "estimate_flow_demo" and other refs in
%                           ijcv_flow code, Sun et al, 2014
%                       Default = 'classic+nl-fastp'
%               inScale = [Tmin Tmax] - min max values for scaling thermal
%                   images to [0 255].
%               flowParams = cell of additional optical flow name-value
%                   pair inputs for 'estimate_flow_interface'
%                      e.g. {'lambda',7,'pyramid_spacing',1.5}
%
%
%          ofile = name PREFIX for output data cube. If empty, velocity
%                   data WILL NOT be saved to discuss, only output into 
%                   workspace by the function
%
% OUTPUT:
%         V = struct containing optical flow results with FIELDS:
%               dataCube = path to therm
%               Vx = x component of velocity field
%               Vz = z (vertical) component of velocity field
%               x  = x spatial vector (meters)
%               z  = z spatial vector (meters)
%               dx = x spacing (meters)
%               dz = z spacing (meters)
%       
% % C Rowell Mar 2020
% Applying Sun et al., (2014) optical flow algorithm to thermal imagery.
% Follwing similar methods to Tournigand et al., (2017), using the 
% ijcv_flow_code package by Sun et al., 2010
%
% 2O: Testing out an approach that approximates a 2nd order central 
% difference scheme...

fprintf('\n========= Thermal Optical Flow Analysis 2nd Order =========\n')
clear textprogressbar
%%
    tic
    disp('Loading data cube...')
    load(thermCube)
    
    % DEFAULTS
    opticDef.Sub = 1:size(D.T,3);
    opticDef.maskPad = 64;
    opticDef.inScale = []; 
    opticDef.method = 'classic+nl-fastp';
    opticDef.flowParams = {};
    
    if nargin<3
        ofile=[];
    end
    if nargin<2
        opticParams = opticDef;
    end
    if isempty(opticParams)
        opticParams = opticDef;
    end
    
    % Check for param fields, set defaults as needed
    fn = fieldnames(opticDef);
    for ff = 1:length(fn)
        if ~isfield(opticParams,fn{ff})
            opticParams.(fn{ff}) = opticDef.(fn{ff});
        elseif isempty(opticParams.(fn{ff}))
            opticParams.(fn{ff}) = opticDef.(fn{ff});
        end
    end
%   Going now by data cube subscript, not frame index
    idx = opticParams.Sub;

        %% Do the thing - setup params and get some stats and initial scaling


    N = numel(idx);
    dii = ismember(D.idx,string(idx)); % Find frame indices

    V.dataCube = thermCube;
    V.opticParams = opticParams;
    V.idx      = idx;
    V.t        = D.t(idx);
    V.ROI      = D.ROI;
    V.dx       = D.dx;
    V.dz       = D.dz;
    V.x        = D.x;
    V.z        = D.z;
    dt = [diff(V.t(1:2)); diff(V.t(2:end))+diff(V.t(1:end-1)); diff(V.t(end-1:end))];
    
    % Vel data cube
    V.Vx = single(zeros([length(D.z) length(D.x) N]));
    V.Vz = V.Vx;
    Vblank = squeeze(V.Vx(:,:,1));
    textprogressbar(sprintf('Calculating %i velocity fields...',length(idx)))
    textprogressbar(0)
    for kk=1:N
        ix = idx(kk);
        if and(kk~=1,kk~=N) % "Centered difference"
            ix1 = idx(kk-1);
            ix2 = idx(kk+1);
        elseif kk==1 % "Forward difference"
            ix1 = idx(kk);
            ix2 = idx(kk+1);
        elseif kk==N
            ix1 = idx(kk-1);
            ix2 = idx(kk);
        end
        [Vroi,~,~] = getROI(any(D.mask(:,:,ix1:ix2),3), 'minSize', 128, 'pad', opticParams.maskPad);

        if ~isempty(Vroi)
            Frame0 = double(D.T(Vroi(1):Vroi(2),Vroi(3):Vroi(4),ix1));
            Frame1 = double(D.T(Vroi(1):Vroi(2),Vroi(3):Vroi(4),ix2));
            
            if ~isempty(opticParams.inScale)
                Frame0 = scaleFrame(Frame0,min(opticParams.inScale),max(opticParams.inScale));
                Frame1 = scaleFrame(Frame1,min(opticParams.inScale),max(opticParams.inScale));
            end
            if ~isempty(opticParams.flowParams)
                uv = estimate_flow_interface(Frame0, Frame1, opticParams.method, [], opticParams.flowParams);
            else
                uv = estimate_flow_interface(Frame0, Frame1, opticParams.method); %, [], varargin);
            end
        
            V.Vx(Vroi(1):Vroi(2),Vroi(3):Vroi(4),kk) = uv(:,:,1);
            V.Vz(Vroi(1):Vroi(2),Vroi(3):Vroi(4),kk) = uv(:,:,2);
        end
        textprogressbar((kk)/(N)*100)
    end
    textprogressbar(' --> Done!')
    
    disp('Converting to real velocities...')
    % Conversion to velocity, flipud to revert to data reference frame
    V.Vx = V.Vx*V.dx./reshape(dt, [1 1 N]);
    V.Vz = V.Vz*V.dz./reshape(dt, [1 1 N]); 
    
    if ~isempty(ofile)
        [vpath,~,~] = fileparts(thermCube);
        vfile = fullfile(vpath, sprintf('%s_%s_n%i.mat', ...
            ofile,datestr(now,'YY-mm-dd'), size(V.Vx,3)));
         fprintf('Optic flow output file:\n\t%s\n',vfile)
        disp('Writing output to disk...')
        save(vfile,'V','opticParams','-v7.3')
    end
    toc
end

function imo = scaleFrame(im,ilow,ihigh,olow,ohigh)
% Same as 'scale_image' in ijcv flow code, but explicitly sets pixels
% outside range to their min/max values

if nargin < 2
    ilow    = min(im(:));
end
if nargin < 3
    ihigh   = max(im(:));
end
if nargin < 4
    olow = 0;
end
if nargin < 5
    ohigh = 255;
end

imo = (im-ilow)/(ihigh-ilow) * (ohigh-olow) + olow;


imo(imo<olow) = olow;
imo(imo>ohigh) = ohigh;


end

