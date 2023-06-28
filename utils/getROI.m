function [roi, omask, poly] = getROI(mask,varargin)
% [roi, omask, poly] = getROI(mask,varargin)
% Generates region of interest bounds from an input plume mask plus height
% limits.
% INPUT:
%   mask  = logical plume mask. Will trim to bounding box plus pad
% OPTIONAL NAME/VALUE PAIRS:
%   index = vector of index locations to overwrite true values in mask. A
%           blank mask with dimensions equal to input mask is created, with
%           values in index vector set to true.
%   iLims = top/bottom limits, pixel units. Can be any length vector of i 
%           subscripts in image. The min and max values will be taken as
%           bounds.
%   jLims = left/right limits, pixel units. Same format as "iLims"
%   pad        = pixel pad around region of interest
%   minSize    = minimum dimensions: scalar or vector of length 2 [i j]
%   maxRegions = scalar integer: allow only N (usually 1) largest contiguous 
%               regions in output mask/poly. Default = []
%   polyMode   = char: 'bounding' or 'separate'; 
%                  -> 'separate' gives separate polygons for each separate
%                      region of true pixels. [Default]
%                  -> 'bounding' gives a single convex polygon bounding all
%                      true pixels
%   crop       = true to crop output mask to roi limits. Default false.
%
%
% OUTPUT:
%   roi  = [i1 i2 j1 j2] aka [y1 y2 x1 x2]
%   mask = mask of same size as original, but cropped to only true values
%           inside roi
%   poly = polygon(s) outlining true pixels
%
% C Rowell April 2020

%% Parsing
    narginchk(1,inf)

    defpad      = 0;
    defIndex    = [];
    defiLims    = [0 0];
    defjLims    = [0 0];
    defMinSz    = 0;
    defMaxReg   = [];
    defPolyMode = 'separate';

    p = inputParser;
    addParameter(p,'pad',defpad)
    addParameter(p,'index',defIndex,@(x)validateattributes(x,{'numeric'},{'vector','integer'}))
    addParameter(p,'iLims',defiLims,@(x)validateattributes(x,{'numeric'},{'vector','integer'}))
    addParameter(p,'jLims',defiLims,@(x)validateattributes(x,{'numeric'},{'vector','integer'}))
    addParameter(p,'minSize',defMinSz,@(x)validateattributes(x,{'numeric'},{'vector','integer'}))
    addParameter(p,'maxRegions',defMaxReg,@(x)validateattributes(x,{'numeric'},{'scalar','integer'}))
    addParameter(p,'polyMode',defPolyMode)
    addParameter(p,'crop',false)
    parse(p,varargin{:})

    opts = p.Results;

%     assert(and(isvector(opts.iLims),length(opt.siLims)>=2),'Input iLims must be vector with length >= 2')
%     iLims = uint8(iLims);

    % Overwrite mask with indices?
    if ~isempty(opts.index)
        mask = false(size(mask));
        mask(opts.index) = true;
    end
%% Get mask bounds

    % Check to make sure there are true values, return if not
    if ~any(mask(:))
        roi = [];
        omask = mask;
        poly = [];
        return
    end

    szy = size(mask,1);
    szx = size(mask,2);
%     yyxx0 = [1 szy 1 szx];
    
    % Initial mask bounding box
    roi(1) = find(sum(mask,2),1,'first');
    roi(2) = find(sum(mask,2),1,'last');
    roi(3) = find(sum(mask,1),1,'first');
    roi(4) = find(sum(mask,1),1,'last');    
%% iLims?
    if all(opts.iLims)
        roi(1) = max([1 roi(1) nanmin(opts.iLims)]);
        roi(2) = min([szy roi(2) nanmax(opts.iLims)]);        
    end
%% jLims?
    if all(opts.jLims)
        roi(3) = max([1 roi(3) nanmin(opts.jLims)]);
        roi(4) = min([szx roi(4) nanmax(opts.jLims)]);        
    end
%% Get output mask
    omask = false(size(mask));
    omask(roi(1):roi(2),roi(3):roi(4)) = mask(roi(1):roi(2),roi(3):roi(4));

    % Check for multiple continguous regions?
    if ~isempty(opts.maxRegions)
        CC = bwconncomp(omask);        
        if length(CC.PixelIdxList)>opts.maxRegions
            [~,biggestI] = sort(cellfun(@numel,CC.PixelIdxList),'descend');
            keepI = biggestI(1:opts.maxRegions);
            
            omask = false(size(mask));
            omask(CC.PixelIdxList{keepI}) = true;
        end
    end
    
    % Update roi
    roi(1) = find(sum(omask,2),1,'first');
    roi(2) = find(sum(omask,2),1,'last');
    roi(3) = find(sum(omask,1),1,'first');
    roi(4) = find(sum(omask,1),1,'last');      
%% Pads?
    if opts.pad>0
        roi(1) = max([roi(1)-opts.pad 1  ]);
        roi(2) = min([roi(2)+opts.pad szy]);
        roi(3) = max([roi(3)-opts.pad 1  ]);
        roi(4) = min([roi(4)+opts.pad szx]);
    end

%% Check Size
    if ~isempty(opts.minSize)
        % Check for minimumm size
        ydim = (roi(2)-roi(1)+1);
        xdim = (roi(4)-roi(3)+1);
        
        if numel(opts.minSize)==1
            minSzX = opts.minSize;
            minSzY = minSzX;
        else
            minSzX = opts.minSize(2);
            minSzY = opts.minSize(1);
        end
        
        if ydim<minSzY
            miny = max([1 round(roi(1)-(minSzY/2-ydim/2))]);
            maxy = min([szy (miny + minSzY - 1)]);
            roi(1:2) = [miny maxy];
        end
        if xdim<minSzX
            minx = max([1 round(roi(3)-(minSzX/2-xdim/2))]);
            maxx = min([szx (minx + minSzX - 1)]);
            roi(3:4) = [minx maxx];
        end
    end
%% Get output polygon bounding actual true values
    switch opts.polyMode
        case 'separate'
            poly = mask2poly(omask);
            
            % Check # Regions
            if ~isempty(opts.maxRegions)
                if length(poly)>opts.maxRegions
                    [~,pI] = sort([poly.Length],'descend');
                    poly = poly(pI(1:opts.maxRegions));
                end
            end
        case 'bounding'
            ch = bwconvhull(omask);
            poly = mask2poly(ch);
    end

    if opts.crop
        omask = omask(roi(1):roi(2),roi(3):roi(4));
    end
    
end