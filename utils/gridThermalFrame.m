function [Frame,mask,gx,gz,dx,dz] = gridThermalFrame(Frame,mask,x,z,dxdz,roundVal)
%   [Frame,mask,gx,gz,dx,dz]
%   Use gridded interpolant to remap a thermal image onto regular grid
%   after mapping pixels to meters.
%
%   IN:     Frame
%           mask
%           x: output of pixelx > meshgrid > px2m
%           z: same as x
%           dxdz = [ds] (equal pixel spacing), [dx dz] (unequal pixel spacing)
%               --> Enter scalar 0 here to force dz=dx, automatically determining
%                   the value
%               --> [0 0] allows dx and dz to be calculated separately
%               --> Non-zero scalar to force both dx,dz to this value
%               --> Non-zero [1 x 2] to force different values
%           roundVal = round dx,dz to nearest roundVal (default: assumes
%           units meters, roundVal = 0.1) 
%
%   OUT:    gFrame
%           gMask
%           gx
%           gz
%           dx
%           dy
%
%   C Rowell Dec 2018

if nargin<6
    roundVal = [];
end
if nargin<5
    dxdz = [];
end
if isempty(roundVal)
    roundVal = 0.1; % Round to nearest? Trying that for now...
end
if isempty(dxdz)
    dxdz = 0;
end

% ---- Default calculated values ------
dx0 = diff(x,1,2); dx0 = round(median(dx0(:))/roundVal)*roundVal;
dz0 = diff(z,1,1); dz0 = round(median(dz0(:))/roundVal)*roundVal;
sign_dx = sign(dx0);
sign_dz = sign(dz0);

% vv % Sets dx and dz equal. Going for "round" of mean to try to mediate between
% losing information in the higher rez axis (likely x) and over-sampling in
% lower rez
ds  = round(mean(abs([dx0 dz0]))./roundVal).*roundVal;

%------ Parse dx,dz values ------

if and(isscalar(dxdz),dxdz==0)
    dx = ds.*sign_dx;
    dz = ds.*sign_dz;
elseif and(isscalar(dxdz),dxdz~=0)
    dx = dxdz;
    dz = dxdz;
elseif numel(dxdz)==2
    if dxdz(1)==0
        dx = dx0;
    else
        dx = dxdz(1);
    end
    if dxdz(2)==0
        dz = dz0;
    else
        dx = dxdz(1);
    end    
else
    error('Input "dxdz" must be scalar or 2-element vector.')
end
if sign(dx)~=sign_dx; dx = -dx; end
if sign(dz)~=sign_dz; dz = -dz; end

%

%%
% get extents for interpolation and new pixel sizes
x0 = max(min(x,[],2));
x1 = min(max(x,[],2));
z0 = min(z(:));
z1 = max(z(:));
% To get new pixel size, interpolate to median value WITHIN mask to
% minimize distortion


gx = x0:dx:x1; gz = z1:dz:z0;
[xq,zq] = meshgrid(gx,gz);

%>>>> RUN INTERPOLANT
FI = scatteredInterpolant(x(:),z(:),Frame(:)); 
Frame = round(FI(xq,zq),2); % Round to 2 decimals to reduce storage size
% mask = MI(xq,zq);
if ~isempty(mask)
    MI = scatteredInterpolant(x(:),z(:),double(mask(:)),'nearest');
    mask  = logical(MI(xq,zq));
end


end
