function [z,R,satfrac,filtfrac] = getRcurve(tk,Rcurve,idx,zLims,zCurve) %,sLength)
% tk = track struct
% Rcurve = choice of radius measure
%               'npx'           Cluster area effective radius, single frame
%               'npxS'            '-> Averaged over window length
%               'npxInt'          '-> Single frame, Interpolated to image z
%               'npxSInt'         '-> windowed average, Interpolated to image z
%               'Tgauss'        Plume Temperature image gaussian coefficients
%               'Vgauss'        Plume Velocity image gaussian coefficients
%               'TkGauss'       Track image gaussian coefficients
%               'plumeWidth'    Width of plume image mask
%               'trackWidth'    Width of track image mask
%               'combined'      Averages npxInt, Tgauss, Vgauss measures
%               'trkCombined'   Averages npxInt and trackWidth measures
% idx    = indices for subsetting R and z curves. [2 x 1] limits or vector
%           of specific indices
% zLims  = lower/upper limits [2 x 1] for z curve. Overrides idx
% Zcurve = choice of z curve for cluster trackStats, 'z' (mean height) or
%           'zmax' (max height). What about a T-weighted mean?
% sLength = smooth length?


if nargin<5 || isempty(zCurve)
    zCurve = 'z';
end
if nargin<4
    zLims = [];
end
if nargin<3
    idx = [];
end
if nargin<2 || isempty(Rcurve)
    Rcurve = 'npx'; % Default for now
end

assert( or(strcmp(zCurve,'z'), strcmp(zCurve,'zmax')),'zCurve not recognized' )
if ~contains(Rcurve,'npx')
    zCurve = 'z';
end


    % Get curves
    switch Rcurve
        case 'npx' 
            R = tk.dat.trackStats.rNpx; 
            %sqrt(tk(kk).npx./pi).*mean(diff(tk(kk).dat.plumeImg.z));
            z = tk.dat.trackStats.(zCurve);
            filtfrac = tk.dat.trackStats.filtFraction;
            satfrac  = tk.dat.trackStats.saturation./tk.dat.trackStats.N;
            
        case 'npxS'
            R = tk.dat.trackStatsS.rNpx; 
            z = tk.dat.trackStatsS.(zCurve);
            filtfrac = tk.dat.trackStatsS.filtFraction;
            satfrac  = tk.dat.trackStatsS.saturation./tk.dat.trackStatsS.N;
             
        case 'npxInt'
            R = tk.dat.trackStatsInt.rNpx; 
            z = tk.dat.trackStatsInt.(zCurve);
            filtfrac = tk.dat.trackStatsInt.filtFraction;
            satfrac  = tk.dat.trackStatsInt.saturation./tk.dat.trackStatsInt.N;
             
        case 'npxSInt'
            R = tk.dat.trackStatsSInt.rNpx; 
            z = tk.dat.trackStatsSInt.(zCurve);
             filtfrac = tk.dat.trackStatsSInt.filtFraction;
            satfrac  = tk.dat.trackStatsSInt.saturation./tk.dat.trackStatsSInt.N;
            
        case 'plumeWidth'
            mi = tk.dat.plumeImg.maskI;
            z = tk.dat.plumeImg.z(mi);
            R = sum(tk.dat.plumeImg.mask,2).*mean(diff(z))/2;
             filtfrac = tk.dat.plumeImg.filtFraction(mi);
            satfrac  = tk.dat.plumeImg.saturation(mi)./tk.dat.plumeImg.N(mi);
           
        case 'trackWidth'
            mi = tk.dat.trackImg.maskI;
            z = tk.dat.trackImg.z(mi);
            R = sum(tk.dat.trackImg.mask(mi,:),2).*mean(diff(z))/2;
            filtfrac = tk.dat.trackImg.filtFraction(mi);
            satfrac  = tk.dat.trackImg.saturation(mi)./tk.dat.trackImg.N(mi);
            
        case 'Tgauss'
            mi = tk.dat.plumeImg.maskI;
            z = tk.dat.plumeImg.z(mi);
            R = tk.dat.plumeImg.gaussCoeffs(:,3);
            filtfrac = tk.dat.plumeImg.filtFraction;
            satfrac  = tk.dat.plumeImg.saturation(mi)./tk.dat.plumeImg.N(mi);
            
        case 'TkGauss'
            mi = tk.dat.trackImg.maskI;
            z = tk.dat.trackImg.z(mi);
            R = tk.dat.trackImg.gaussCoeffs(:,3);
            filtfrac = tk.dat.trackImg.filtFraction(mi);
            satfrac  = tk.dat.trackImg.saturation(mi)./tk.dat.trackImg.N(mi);
            
        case 'Vgauss'
            mi = tk.dat.plumeVz.maskI;
            z = tk.dat.plumeVz.z(mi);
            R = tk.dat.plumeVz.gaussCoeffs(:,3);
             filtfrac = tk.dat.plumeVz.filtFraction(mi);
            satfrac  = tk.dat.plumeImg.saturation(mi)./tk.dat.plumeImg.N(mi);
           
        case 'combined'
            mi = tk.dat.plumeImg.maskI;
            z1 = tk.dat.plumeImg.z(mi);
            z2 = tk.dat.trackStatsInt.z; 
            
            [z,ia,ib] = intersect(z1,z2);
            R = zeros(length(z),3); 
            sf = R;
            ff = R;
            [~,R(:,1),sf(:,1),ff(:,1)] = getRcurve(tk,'Tgauss',ia);
            [~,R(:,2),sf(:,2),ff(:,2)] = getRcurve(tk,'Vgauss',ia);
            [~,R(:,3),sf(:,3),ff(:,3)] = getRcurve(tk,'npxInt',ib); % Optionally smooth this curve?
            R        = mean(R,2);
            filtfrac = mean(ff,2);
            satfrac  = mean(sf,2);
            
        case 'trkCombined' % Using gauss typically does not make sense for tracks
            mi = tk.dat.trackImg.maskI;
            z1 = tk.dat.trackImg.z(mi);
            z2 = tk.dat.trackStatsInt.z; 
            [z,ia,ib] = intersect(z1,z2);
            R = zeros(length(z),2); 
            sf = R;
            ff = R;
            [~,R(:,1),sf(:,1),ff(:,1)] = getRcurve(tk,'trackWidth',ia);
            [~,R(:,2),sf(:,2),ff(:,2)] = getRcurve(tk,'npxInt',ib);
            R        = mean(R,2);
            filtfrac = mean(ff,2);
            satfrac  = mean(sf,2);
            
        otherwise
            error('R curve not recognized')
    end

    % Subsetting curves
    if ~isempty(zLims)
        if isnan(zLims(1)); zLims(1) = min(z); end
        if isnan(zLims(2)); zLims(2) = max(z); end
        assert(numel(zLims)==2,'zLims must by a vector of length(2)')
        zi          = and( z>min(zLims), z<max(zLims) );
        z           = z(zi);
        R           = R(zi);
        filtfrac    = filtfrac(zi);
        satfrac     = satfrac(zi);
        
    elseif ~isempty(idx)
        if numel(idx)==2
            if isnan(idx(1)); idx(1) = 1;           end
            if isnan(idx(2)); idx(2) = length(z);   end
            z        = z(min(idx):max(idx));
            R        = R(min(idx):max(idx));
            satfrac  = satfrac(min(idx):max(idx));
            filtfrac = filtfrac(min(idx):max(idx));
            
        else
            assert(length(idx) <= length(z),'length(idx) must be <= length(z) or length(2)')
            z        = z(idx);
            R        = R(idx);
            satfrac  = satfrac(idx);
            filtfrac = filtfrac(idx);
        end
        
    end
    
end