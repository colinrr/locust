function [zErr,zErrMax,Rstat] = zErrorEstimation(x,z,mask,geom,metric)
% [zErr,zErrMax] = zErrorEstimation(x,z,mask,geom,metric)
% ESTIMATE vertical position error in plume imagery, assuming a vertical
% projection plane.
%   x       = x position vector or matrix matching T RELATIVE TO CAMERA CENTER PIX
%            (will be centered to vent internally)
%   z       = z position vector or matrix matching T RELATIVE TO CAMERA ELEVATION
%   T       = 2- or 3D array of shape masks
%   geom    = mapPixels output geometry structure
%   metric  = 'cylindrical' or 'bwdist' (Default = cylindrical). Bwdist
%           should only be used on a regular gridded image. If the x and z
%           coordinates are not regular, use cylindrical.
%
% OUT:
%   zErr    = image of estimated z error based on plume radius
%   zErrMax = Range of z error uncertainty based on plume axial position
%   Rstat   = error geometry parameters
%
% C Rowell Jun 2020

    if nargin<5
        metric = 'cylindrical';
    end
    
    N = size(mask,3);

    assert(or(all(size(x)==size(mask,[1 2])), numel(x)==size(mask,2)),'Size of x does not match mask')
    assert(or(all(size(z)==size(mask,[1 2])), numel(z)==size(mask,1)),'Size of z does not match mask')
    assert(any(size(x,3)==[N 1]),'3rd dimension of x must be 1 or size(mask,3)')
    assert(any(size(z,3)==[N 1]),'3rd dimension of z must be 1 or size(mask,3)')

    zErr = zeros(size(mask));
    zErrMax = zeros(size(mask,1),3,N);
    
    xvec = and(isvector(x),numel(x)==size(mask,2));
    zvec = and(isvector(z),numel(z)==size(mask,1));
    if and(xvec,zvec)
        [x,z] = meshgrid(x,z);
    elseif xvec
        x = repmat(reshape(x,[1 length(x)]),[size(mask,1) 1]);
    elseif zvec
        z = repmat(reshape(z,[length(z) 1]),[size(mask,2) 1]);
    end

    [xV,zV] = px2m(geom.target_pix_ij(2),geom.target_pix_ij(1),geom);
    x = x - xV;  % Center x over vent
    theta = atand(z./geom.X_distance); % Elevation angle at each point
    
    for nn=1:N
        maskn = mask(:,:,nn);
        if any(maskn(:))
            switch metric
                case 'cylindrical'
                    try maskn = trimMask(maskn,2);
                    catch me
                        warning('trimMask failure in zErrorEstimation, frame %i',nn)
                    end
            end
            
            % Get trimmed mask bounds
            npix = sum(maskn,2);
            pfirst = find(npix>0,1,'first'); % plume MASK top row
            plast = find(npix>0 ,1,'last');  % plume MASK bottom row
            plvec = find(any(maskn,2))';
            masklims = zeros(size(maskn,1),2);
            for mm = plvec %pfirst:plast
                masklims(mm,1) = find(maskn(mm,:),1,'first');
                masklims(mm,2) = find(maskn(mm,:),1,'last');
            end
            noMask = any(~masklims,2);
            nsmooth = max([7 round(median(masklims(~noMask,2)-masklims(~noMask,1))/mean(diff(x(1,:)))/4)]);
            masklims(:,1) = round(smooth(masklims(:,1),nsmooth));
            masklims(:,2) = round(smooth(masklims(:,2),nsmooth));
            masklims(noMask,:) = NaN;

%             leftI   = sub2ind(size(maskn),(pfirst:plast)',masklims(~noMask,1));
%             rightI  = sub2ind(size(maskn),(pfirst:plast)',masklims(~noMask,2));
            try leftI   = sub2ind(size(maskn),plvec',masklims(~noMask,1));
            catch wtf
                nn
            end
            rightI  = sub2ind(size(maskn),plvec',masklims(~noMask,2));

            R            = zeros(size(maskn,1),1);
            R(~noMask)   = (x(rightI) - x(leftI))/2;
            x_c          = zeros(size(R));
            z_c          = zeros(size(R));
            x_c(~noMask) = x(leftI) + R(~noMask);
            x_c(noMask)  = NaN;
    %         z_c(~noMask) = z(
            dxR = x-x_c;

            switch metric
                case 'bwdist'
                    % Find top and bottom cut limits for mask
                    dz = mean(diff(z(:,1))); % Assumes regular z!
                    del2R = del2(smooth(R,3),dz);


                    cutlength = round(sum(~noMask)*0.05);    
                    if pfirst>1
                        [~,Rcut1] = min(del2R(pfirst:pfirst+cutlength));
                        Rcut1 = Rcut1 + pfirst - 1;
                    else
                        Rcut1 = 1;
                    end
                    if plast<size(maskn,1)
                        [~,Rcut2] = min(del2R(plast-cutlength:plast));
                        Rcut2 = Rcut2 + plast-cutlength-1;
                    else
                        Rcut2 = size(maskn,1);
                    end

                    bwmask = zeros(size(mask(:,:,nn)));
                    bwmask(leftI) = 1; bwmask(rightI) = 1;
                    bd = bwdist(~mask(Rcut1:Rcut2,:,nn));
                    [Rbd,Rj] = max(bd,[],2);
                    R(Rcut1:Rcut2) = Rbd*abs(dz);
                    Rix = sub2ind(size(maskn),(1:size(bd,1))'+pfirst-1,Rj);
                    x_c(Rcut1:Rcut2) = x(Rix);

                    dxR(Rcut1:Rcut2,:) = (Rbd-bd)*abs(dz);
            end

        % y_c must be set or assumed. Opts 2,3 give a decent estimate of full range
        y_c = 0;
    %     y_c = abs(x_c); % Option 2 
    %     y_c = -abs(x_c); % Option 3

            dz_a = tand(theta) .* (sqrt(R.^2-(dxR).^2) + y_c);
            zErr(:,:,nn) = real(dz_a).*maskn;

    %         dZmax = zeros(size(dz_a,1),3);
            zErrMax(:,1,nn) = tand(theta(:,1)) .* (R - abs(x_c));
            zErrMax(:,2,nn) = tand(theta(:,1)) .* (R);
            zErrMax(:,3,nn) = tand(theta(:,1)) .* (R + abs(x_c));
        end
    end

    Rstat = [];
    % Temporary output
    if N == 1
        Rstat.xV = xV;
        Rstat.R = R;
        Rstat.xc = x_c+xV;
        Rstat.theta = theta;
    else
        Rstat = [];
    end

    
end