function [x,z] = px2m(jj,ii,g,projection)
% [x,z] = px2m(px,pz,g,projection)
% IN:   jj = input x (horizontal) vector or matrix (pixel coords, dimension 2)
%       ii = input z (vertical) vector or matrix (pixel coords, dimension 1)
%       g  = camera projection geometry structure. Required fields:
%               X_distance
%               nx0
%               ny0
%               Elev_angle
%               hifov_px_deg
%               vifov_px_deg
%       proejction = 'rectilinear' or 'original' (for legacy purposes)
%
% OUT:  x = output x vector or matrix, meters
%       z = output z vector or matrix, meters

if nargin<4
    projection = 'rectilinear';
end

% Check for registration and adjust pixel coordinates
if isfield(g,'unregistered')
%     fprintf('NOTE: Correcting x,z pixels for registered images.\n')
    jj = jj + g.regXlim(1) - 1;
    ii = ii + g.regYlim(1) - 1;
end

switch projection
    case 'rectilinear' % Using new geometry struct parameters
        z = g.X_distance*( tand(g.Elev_angle - ((ii-0.5)-g.ny0)*g.vifov_px_deg) );
        x = (g.D_los + sind(g.Elev_angle).*(z-g.X_distance.*tand(g.Elev_angle))).*tand(((jj-0.5)-g.nx0)*g.hifov_px_deg);
        
        
        
%     case 'rectilinear_legacy' % Using old geometry struct parameters
        
    case 'original' % Initial projection used in 2018.
        x = g.X_distance* tand(((jj-0.5)-g.nx0)*g.hifov_px_deg)./cosd(g.Elev_angle - (ii-0.5-g.ny0)*g.vifov_px_deg);
        z = g.X_distance*( tand(g.Elev_angle - ((ii-0.5)-g.ny0)*g.vifov_px_deg) );
    otherwise
        error('Projection type for pixel mapping not recognized')
end

end