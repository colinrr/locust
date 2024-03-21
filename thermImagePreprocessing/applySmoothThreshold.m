function maskThermal(inputDir,matHeads,Idx,outputDir,ITW,mask,nullVal)
% maskThermal(inputDir,matHeads,Idx,outputDir,ITW,mask,nullVal)
% Applies smooth-heaviside thresholding filter and/or mask to frames,
% outputs modified frames into new directory. Intended for pre-processing a
% copy of thermal imagery prior to use in Bombrun et al. (2018) "plumeTracker".
%   Usefull tools to help design this:
%       makeForegroundMask
%       getTimeHistogram
%       plotTimeHistogram
%
%   inputDir = directory of .mat files with image data
%   heads    = header table file associated with image data
%
%   Idx     = list of frame indices to apply to. Defaults to all in "heads"
%   outputDir  = [default: inputDir/preProd/ ]
%   ITW =   Secondary threshold to kill background values. Must be [3 x n].
%           Because objects of interest or background temperature values may change 
%           with time, this allow a time-varying threshold function to be applied to
%           thermal sequence. The filter will be interpolated between these values.
%           Use getTimeHistogram + plotTimeHistogram to visualize
%       eg.  ITW = [1 1000 2400 ; ... % Frame Index for control points.
%                   255 245 240 ; ... % T0 - Center of smooth threshold filter, Kelvin
%                   10  10  10  ];    % Width of smooth threshold filter, in Kelvin
%   mask   = logical mask, dimensions = size(Frames). 
%           Frame pixels with true mask values will be set to nullVal
%   
% could implement maskIdx to select frames for applying masks if needed.
% C Rowell March 2020

if nargin<7
    nullVal = [];
end
if nargin<6
    maskFile = '';
end
if nargin<5
    ITW = [];
end
if and(isempty(ITW),isempty(maskFile))
    error('Needs input of either a mask or temperature threshold')
end

%% Setup
load(matHeads)
if isempty(Idx)
    Idx = str2double(T.Properties.RowNames); % String or double?
end

if ~isempty(ITW)
    apply_heaviside = true;
    %Threshold function
    hs_sm = (@(T,T0,k) 1./(1 + exp(-2*k*(T-T0))));
    
    % Get threshold vectors
    if size(ITW,2)>1
        T0v = interp1(ITW(1,:),ITW(2,:),Idx,'pchip');
        dTv = interp1(ITW(1,:),ITW(3,:),Idx,'pchip');
    else
        T0v = ITW(2)*ones(size(Idx));
        dTv = ITW(3)*ones(size(Idx));
    end
    Kv = 2./dTv;
else
    apply_heaviside = false;
end

if ~isempty(maskFile)
    apply_mask = true;
    load(maskFile)
    assert(islogical(mask),'input "mask" must be logical.')
end
%% Loop over frames
fprintf('Editing %i input frames..\n.',numel(Idx))
if apply_heaviside; disp('Applying threshold filter'); end
if apply_mask; disp('Applying mask(s)'); end
fprintf('Writing files to:\n\t%s\n',outputDir)
dispstat('','init')
for kk = 1:length(Idx)
    
    idx = Idx(kk);
    load(fullfile(inputDir,T.File{num2str(idx)}))
    if and(isempty(nullVal),kk==1)
        nullVal = min(Frame(:));
    end
    if apply_mask
        Frame(mask) = nullVal;
    end
    if apply_heaviside
        Frame = hs_sm(Frame,T0,K) .* (Frame - nullVal) + nullVal;
    end
    oname = fullfile(outputDir,T.File{num2str(idx)});
    dispstat(sprintf('%i/%i: %s',kk,numel(Idx),T.File{num2str(idx)}))
    save(oname,'Frame','File_DateTime')
end

    disp('Writing pre-processing paramter file.')
    preProc.inputDir    = inputDir;
    preProc.Headers     = matHeads;
%     preProc.satVal      = satVal;
    preProc.nullVal     = nullVal;
%     preProc.thresh_fg   = thresh_fg;
%     preProc.ref_idx     = ref_idx;
    preProc.Idx         = Idx;
%     preProc.peakRange   = peakRange;
    preProc.ITW         = ITW;
%     preProc.interp_meth = interp_meth;
%     preProc.bins        = bins;
    preProc.mask        = mask;
    if ~isempty(ITW)
        preProc.T0v = T0v;
        preProc.dTv = dTv;
        preProc.Kv  = Kv;
    end
%     save(fullfile(outputDir,'preProcParams.mat'),'preProc');
end