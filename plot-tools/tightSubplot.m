function varargout = tightSubplot(nrows,ncols,N,dx,dy,pads,xsize,ysize)
% [ax] = tightSubplot(nrows,ncols,N,dx,dy,pads,xsize,ysize);
% Slightly customized subplotting routine for cleaner/tighter axes than a
% standard subplot.
% INPUT:    nrows = same as in subplot
%           ncols = same as in subplot
%           N     = same as in subplot
%         Optional:
%           dx    = normalized x spacing between axes, 0 to 1 [adaptive default
%                   - 25% current axis size] [ACTUALLY just defaults to 0 now]
%           dy    = normalized y spacing between axes, 0 to 1 [adaptive default
%                   - 30% current axis size]
%           pads  = normalized border space vector [left right bottom top]
%                   Default: [0.11 0.04 0.1 0.06]
%           xsize = vector that generates non-uniform axes widths,
%                   normalized to each other
%           ysize = vector that generates non-uniform axes heights,
%                   normalized to each other
%      
%           Setting any optional argument to 0 will turn off
%           ticks/labels on the corresponding axis. 
%           Eg 1a: dx = 0 will turn off axis labels and tick labels for the y
%           axis, EXCEPT for the left most axis, since the plots will have
%           no spacing between them.
%           Eg 1b: pads = [0 0.1 0.1 0.1]  would ALSO turn off the left most
%           axis, since you are leaving no border space on the left side of
%           the figure.
%   
%           Eg 2. xsize = [2 1] will tile columns of subplots with RELATIVE
%           width ratios of 2 to 1, filling space allotted by pads and dx.
%           
%               tightSubplot(4,2,N,[],[],[],[2 1]) for N = 1:8 will create
%           a plot that looks like this:
%                   
%                       |     |   |     |
%                       |     |   |     |
%                       |---- |-- |---- |--
%
%                       |     |   |     |
%                       |     |   |     |
%                       |---- |-- |---- |--
%
%   OUTPUT:  ax = axes handle
%
% C Rowell, Jun 2017
%     NOTES:  Might implement variable dx,dy spacing, adaptive pads based
%     on aspect ratio
%% Set defaults
% MATLAB DEFAULTS
% pads_def = [0.13 0.095 0.11 0.075]; % Matlab defaults
% dx_X_ratio = 0.3158; % Ratio of horizontal subplot spacing to axes width
% dx_X_ratio = 0.3889; % Ratio of vertical subplot spacing to axes height
nargoutchk(0,1)

% My defaults
pads_def = [0.11 0.04 0.1 0.06];
dx_X_ratio = 0; %0.25; % Changed these to 0 for now, makes more sense as the default here
dx_Y_ratio = 0; %0.3;
% Aspect ratio thresholds - experimental
aspR_hi = 1.25; % If higher, use xpads as ypads
aspR_lo = 0.5;  % If lower, use ypads as xpads
%% Parse input
if nargin<4
    dx = [];
end
if nargin<5
    dy = [];
end
if nargin<6
    pads = [];
end
if nargin<7
    xsize = [];
end
if nargin<8
    ysize = [];
end

%% Get aspect ratio, calculate spacings
% This will initiate a figure if one does not exist
figpos = get(gcf,'position');
aspRatio = figpos(4)/figpos(3);

% Assign pads if necessary
if isempty(pads)
    pads = pads_def;
end
% Check aspect ratio
% if aspRatio>=aspR_hi
%     pads(3:4) = pads(3:4)/aspRatio;
% elseif aspRatio<=aspR_lo
%     pads(1:2) = pads(1:2)*aspRatio;
% end
%% Calc spacings where necessary
% Number of spaces needed
ndx = ncols-1;
ndy = nrows-1;

% Calc total space available, and average space used by axes
x_space = 1-sum(pads(1:2));
y_space = 1-sum(pads(3:4));

if isempty(dx)
    dx = (x_space)/ncols*dx_X_ratio; % Calc dx assuming evenly sized axes
%     if aspRatio<=aspR_lo % Adjust for aspect ratio
%         dx        = dx*aspRatio;
%     end
end
if isempty(dy)
    dy = (y_space)/nrows*dx_Y_ratio; % Calc dy assuming evenly sized axes
%     if aspRatio>=aspR_hi % Adjust for aspect ratio
%         dy        = dy/aspRatio;
%     end
end

%% Calc axis size and position based on pads,dx,dy, and size vectors
% Assign current axes position from index
iy = ceil(N/ncols);
ix = N - (iy-1)*ncols;

% Create normalized length vector for axes
if ~isempty(xsize)
    xsize = repmat(xsize,[1 ceil(ncols/numel(xsize))]);
    xsize = xsize(1:ncols)/sum(xsize(1:ncols))*(x_space-ndx*dx);
else
    xsize = ones(1,ncols)*((x_space - (ndx)*dx)/ncols);
end

if ~isempty(ysize)
    ysize = repmat(ysize,[1 ceil(nrows/numel(ysize))]);
    ysize = ysize(1:nrows)/sum(ysize(1:nrows))*(y_space-ndy*dy);
else
    ysize = ones(1,nrows)*((y_space - (nrows-1)*dy)/nrows);
end

% Calc position vector
xw = xsize(ix); %(x_space - (ncols-1)*dx)/ncols; % Axes width
yh = ysize(iy); %(y_space - (nrows-1)*dy)/nrows; % Axes height
x = pads(1) + (ix-1)*dx + sum(xsize(1:ix-1));
y = pads(3)+y_space-(iy-1)*dy-sum(ysize(1:iy));   % pads(3) + (iy-1)*(yh+dy);
pos = [x y xw yh];

%% Add axes
ax=axes('position',pos);
% Text to check position as a useful trick
% text(0.5,0.5,sprintf('%i, %i',ix,iy))
ax_stat = get(gca);

%% Turn off/adjust axes labels when spacing is 0
if dx==0
    box on
    if ix>1
        set(gca,'YTickLabel',[])
    end
%     xt  = get(gca,'XTick');
%     xtl = get(gca,'XTickLabel');
%     xl  = xlim;
    if ix~=ncols && ax_stat.XTick(end)==ax_stat.XLim(2);
        set(gca,'XTick',ax_stat.XTick(1:end-1))
    end
end
if dy==0
    box on
    if iy<nrows
        set(gca,'XTickLabel',[])
    end
    if iy>1
        set(gca,'YTick',ax_stat.YTick(1:end-1))
    end
end

if nargout==1
    varargout{1} = ax;
end
