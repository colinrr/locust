function plotZerrorFrame(x,z,t,mask,zErr,zErrMax,geom,idx,zRef)
% D = thermal cube
% zErr, zErrMax - outputs 1 and 2 of zErrorEstimation
% geom = geometry struct
% idx  = vector of frame subscripts
% zRef = 'camera' or 'vent'

if nargin<6
    zRef = 'camera';
end
if nargin<5
    idx = [];
end
if isempty(idx)
    idx = 1:length(t);
end
if nargin<4
    geom = [];
end
if isempty(geom)
    geom = geom;
end

% Apparent height calculation
[xV,zV] = px2m(geom.target_pix_ij(2),geom.target_pix_ij(1),geom);
dz0 = geom.Ztarg-geom.Z0;

switch zRef
    case 'camera'
        z = z + dz0;
        zl = 'Z [m above camera]';
    case 'vent'
        z = z;
        zl = 'Z [m above vent]';
    otherwise
        error('Height reference string not recognized')
end

if size(zErr,3)==size(mask,3) && length(idx)~=size(mask,3)
    zErr = zErr(:,:,idx);
end
if size(zErrMax,3)==size(mask,3) && length(idx)~=size(mask,3)
    zErrMax = zErrMax(:,:,idx);
end

for nn=1:numel(idx)
    ix = idx(nn);
    P = mask2poly(mask(:,:,ix));
    [~,pi] = max([P.Length]);
    
    figure
    ax(1) = tightSubplot(1,2,1,0.1,[],[],[3 1]);
    surf(x-xV,z,zErr(:,:,nn))
    shading flat
    view([0 0 1])
    colorbar
    colormap(jet(256))
    axis tight equal
    xlabel('X [m from vent]')
    ylabel(zl)
    title(sprintf('Estimated Z error, t = %s, y_c = 0',t(ix)))
    hold on
    plot(x(P(pi).X)-xV,z(P(pi).Y),'m')
    contour3(x-xV,z,zErr(:,:,nn),[0 50 100 150 200 300 400],'w')
    set(gca,'FontSize',12)
    caxis([0 300])
    
    ax(2) = tightSubplot(1,2,2,0.1,[],[],[3 1]);
    plotLineError(zErrMax(:,:,nn),z)
    set(gca,'FontSize',12)
    xlabel('Estimated range of Z error [m]')
    ylabel(zl)
    grid on
    linkaxes(ax,'y')
    axis tight
end

end