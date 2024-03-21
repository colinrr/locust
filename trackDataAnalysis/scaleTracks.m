function [trackFits,fitArrays,opts] = scaleTracks(tk,varargin) %,trackSet,Tcurve,rIndices,Tindices,sensitivityFlag)
% [tk,fitArrays] = scaleTracks(tk,varargin)
% Get scaled track dataset
% tk        = structure of track images output by createTrackImgs, or path to
%               .mat struct file
%
% NAME/VALUE PAIRS:
% trackSet  = optional indices to subset tk (ie which tracks to use?).
%               Defaults to ALL (1:length(tk)).
% Rcurve    = switch to determine which Radius curve to use for deriving dr/dz
%               and z0. OPTIONS:
%               'npx'           Cluster area effective radius, single frame
%               'npxS'            '-> Averaged over window length
%               'npxInt'          '-> Single frame, Interpolated to image z
%               'npxSInt'         '-> windowed average, Interpolated to image z
%               'Tgauss'        Plume Temperature image gaussian coefficients
%               'Vgauss'        Plume Velocity image gaussian coefficients
%               'TkGauss'       Track image gaussian coefficients
%               'plumeWidth'    Width of plume image mask
%               'trackWidth'    Width of track image mask
%               'combined'      Averages npxInt, Tgauss, Vgauss'npx', 'npxInt' 'Tgauss' ,'Vgauss', 'TkGauss' , 'plumeWidth','trackWidth', 'combined'
%
% Tstruct    'structName' to determine which temperature curve to use
%                 options:          'plumeImg'
%                                   'trackImg'
%                                   'trackStats'
%                                   'trackStatsS'
%                                   'trackStatsInt'
%                                   'trackStatsIntS'
%  
% Tcurve     'curveName' opts:      'prctile' - 95th percentile
%                                   'mean' 
%                                   'max'
%
% Ridx      Optional indices to subset tracks in tk. (ie which parts of
%             each RADIUS track to use?)
%            [N x 2] array (N = length((tk) or (trackSet)): [Ridx1 RidxN] for each track
%               OR
%            [N x 1] cell array of index vectors
%            Use 0 or NaN to default a track to it's normal limits.
%
% Tidx      SAME AS Ridx, applied to T tracks
%
% RzLim     [N x 2] array of min/max z bounds for radius fitting.
% TzLim     [N x 2] array of min/max z bounds for temperature fitting.
%
%           NOTE: indices and z limits are additive. E.g., Specified Ridx 
%           values will be include along with those inside the limits 
%           imposed by RzLim. (?same is true for Tidx/TzLim?)
%
% z0        [N x 3] array of z0 with error estimate. [z0_lowbound z0 z0_hibound]
%               n can be 1, length(tk), or length(trackSet), and will be
%               distributed accordingly to each track.
%             -> Can also be [n x 1], in which case no error bars are
%             calculated.
%             -> Nans in the central estimate (ie second column for [nx3] case)
%             will ignore the manual z0 and revert to one calculated from R values.
%
% r0        [N x 1] vector:  manual normalization value for track heights. 
%
% Bmethod   Method for retrieving temperature power law exponent.
%           'fitPL' (estimate z0 with lin fit - default)
%               OR
%           'George1977' (combine Z0,B estimated) - NOT YET IMPLEMENTED
%
% sensitivityFlag   Retrieves B(z0), the sensitivity curve of the power
%                   law estimate to the choice of z0.
%
% filtFracThresh    Threshold of fraction of pixels filtered above which to ignore
%                   data points for fitting. (default = 0.1, ie points
%                   for which more than 10% of pixels were filtered out
%                   due to height uncertainty will not be included in
%                   the R or T fits).
%
% satFracThresh     Threshold fraction of saturated pixels above which to
%                   ignore data points for fitting. (default = 0.1)
%
% dT0Limits         [Low High] ABSOLUTE bounds (Kelvin) for c coefficient in power law fit
%                   (a*z^b + c). This amounts to max uncertainty in the
%                   background value of deltaT about 0 after atmospheric 
%                   profile removal. Defaults (generously) to +/-0.5 Tmax
%
% cLimits           [Low High] RELATIVE bounds for c coefficient in power
%                   law fit, as a fraction of max(Delta T). The power law
%                   fit will apply the tighter of the absolute and relative
%                   bounds (ie closest to 0).
%
% OUTPUT:  
%   trackFits      Struct of detailed output for each track
%   fitArrays      Struct of array results for quick comparison across tracks
%                   (z0, r0, drdz, B, ?fit quality stats?)

fprintf('\n========= Retrieve Track Scales and Exponents =========\n')

%%  Parse input
assert( isstruct(tk) ,'Input "tk" must be a struct')
% assert( isstruct(tk) || exist(tk,'file') ,'Input "tk" must be a struct or path to struct file')

defDesc             = 'Track results'; % Description of scaled set
defTrackSet         = 1:length(tk); % Which tracks to use?
defUseZmax          = false;        % Use maximum height for tracked features instead of mean? (trackStats curves only)
defRcurve           = 'npxInt';     % Data to use for Radius. OPTIONS: 'npx','npxInt, 'Tgauss' ,'Vgauss', 'TkGauss' , 'plumeWidth','trackWidth', 'combined'
% defZR               = 'trackStats'; %
defZ0               = [];           % Manually enter a z0 (m)?
defR0               = [];           % Manually enter an r0
defTstruct          = 'trackStats'; % Which track structure to use for temperature curve?
defTcurve           = 'prctile';    % Which statistical Temperature curve to use for B estimate?
defRidx             = [];           % Array of integers [ntracks x 2] or [n x 1] cell
defTidx             = [];           % Array of integers [ntracks x 2] or [n x 1] cell
defUseRzLim            = true;           % Array of integers [ntracks x 2] corresponding to min/max heights to fit for each R track
defUseTzLim            = true;           % Array of integers [ntracks x 2] corresponding to min/max heights to fit for each T track
defBmethod          = 'fitPL';      % 'fitPL' (estimate z0 with lin fit - default) or 'George1977' (combine Z0,B estimated)
defSensFlag         = false;        % Run B sensitivity test to get B estimate as a function of Z0
defFiltFracThresh   = 0.1;          % Threshold of filter fraction above which to ignore data points
defSatFracThresh    = 0.1;          % Threshold of saturation fraction above which to ignore data points
defdT0Limits        = [];           % Apply absolute limits (Kelvin) to the c coefficient in the power law fit
defcLimits          = [];           % Apply relative limits (fraction of max deltaT) to the c coefficient in the power law fit
% zLims?    specify height limits for R/T curves instead of indices?

p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'desc',defDesc)
addParameter(p,'trackSet',defTrackSet)
addParameter(p,'useZmax',defUseZmax)
addParameter(p,'Rcurve',defRcurve)
addParameter(p,'Tstruct',defTstruct)
addParameter(p,'Tcurve',defTcurve)
addParameter(p,'z0',defZ0, @(x) (isempty(x) || (isnumeric(x) && or(size(x,2)==3,size(x,2)==1))))
addParameter(p,'r0',defR0, @(x) isnumeric(x) && size(x,1)==1)
addParameter(p,'filtFracThresh',defFiltFracThresh, @(x) isscalar(x))
addParameter(p,'satFracThresh',defSatFracThresh, @(x) isscalar(x))
addParameter(p,'sensitivityFlag',defSensFlag)
addParameter(p,'Bmethod',defBmethod)
addParameter(p,'Ridx',defRidx,@(x) isempty(x) || and(isnumeric(x),size(x,2)==2) || iscell(x))
addParameter(p,'Tidx',defTidx,@(x) isempty(x) || and(isnumeric(x),size(x,2)==2) || iscell(x))
% addParameter(p,'RzLim',defRzLim,@(x) isempty(x) || and(isnumeric(x),size(x,2)==2))
% addParameter(p,'TzLim',defTzLim,@(x) isempty(x) || and(isnumeric(x),size(x,2)==2))
addParameter(p,'useRzLim',defUseRzLim)
addParameter(p,'useTzLim',defUseTzLim)
addParameter(p,'dT0Limits',[])
addParameter(p,'cLimits',[])

parse(p,varargin{:})
opts = p.Results;

Ntracks0    = length(tk);

assert(length(opts.trackSet)<=Ntracks0,'length(opts.trackSet) !<= length(tk) !')

% ------- Parse curve indices or zlimits -------
% R index
if ~isempty(opts.Ridx)
    assert( or(size(opts.Ridx,1)==Ntracks0 , size(opts.Ridx,1)==length(opts.trackSet)),'1st dimension of Indices must match length of tk struct or trackSet')
    if size(opts.Ridx,1)==Ntracks0 && size(opts.Ridx,1)~=length(opts.trackSet)
        opts.Ridx = opts.Ridx(opts.trackSet,:);
    end
    useRidx = true;
else
    useRidx = false;
end


% Assign Z limits - RETRIEVE MANUALLY CHOSEN HEIGHT LIMITS FOR CURVE
% FITTING DATA
if opts.useRzLim || opts.useTzLim
    [RzLim,TzLim,Rqc,Tqc] = getTrackZlimits(opts.trackSet,opts.Rcurve,opts.Tstruct,opts.useZmax);
    if opts.useRzLim
        opts.RzLim = RzLim;
    end
    if opts.useTzLim
        opts.TzLim = TzLim;
    end
else
    Rqc = ones(length(opts.trackSet),1);
    Tqc = ones(length(opts.trackSet),1);
end


% T index
if ~isempty(opts.Tidx)
    assert( or(size(opts.Tidx,1)==Ntracks0 , size(opts.Tidx,1)==length(opts.trackSet)),'1st dimension of Indices must match length of tk struct or trackSet')
    if size(opts.Tidx,1)==Ntracks0 && size(opts.Tidx,1)~=length(opts.trackSet)
        opts.Tidx = opts.Tidx(opts.trackSet,:);
    end
    useTidx = true;
else
    useTidx = false;
end


% Manual Z0 parsing
if ~isempty(opts.z0)
    assert( size(opts.z0,1)==1 || size(opts.z0,1) == Ntracks0 || size(opts.z0,1) == length(opts.trackSet), '1st dimension of manual z0 must be one of: 1, length(tk), length(trackSet)' );
    if size(opts.z0,1) == Ntracks0
        opts.z0 = opts.z0(opts.trackSet,:);
    end
end

% Manual r0 parsing
if ~isempty(opts.r0)
    if size(opts.r0,1) == Ntracks0
        opts.r0 = opts.r0(opts.trackSet);
    end
end


tk          = tk(opts.trackSet);
Ntracks     = length(tk);

%% MAIN TRACK LOOP


% (z0, r0, drdz, B, ?fit quality stats?)
fitArrays.tSpan = zeros(Ntracks,2); % Video time span (s)
fitArrays.z0    = zeros(Ntracks,3); % asymptotic virtual origin guess (m)
fitArrays.r0    = zeros(Ntracks,3); % vent radius guess (m)
fitArrays.drdz  = zeros(Ntracks,3); % drdz slope (-)
fitArrays.r_r2  = zeros(Ntracks,1); % R model adjusted r-square
fitArrays.r_rmse = zeros(Ntracks,1); % R model rmse (m)
fitArrays.Bpl   = zeros(Ntracks,3); % Power law fit
fitArrays.B_r2  = zeros(Ntracks,1); % R^2 for powerlaw fit
fitArrays.B_nrmse = zeros(Ntracks,1); % Normalized rms error for powerlaw fit
% fitArrays.Bfs   = zeros(Ntracks,3); % T vs. z^b search fit
% fitArrays.Blm   = zeros(Ntracks,3); % Linear model fit

% kk = 1;
for kk=Ntracks:-1:1

    % ------- (1) get r(z) -------- 
    

    
    % Get full R curve
    [zR,R,sf,ff] = getRcurve(tk(kk),opts.Rcurve); %,[],[],'z');
    
    % Parse R indices/z limits and apply filter thresholds
    ixfilt = false(length(zR),2);
    sffilt =  ff<=opts.filtFracThresh; % Don't need saturation filtering for r values
    if ~useRidx && ~opts.useRzLim
        ixfilt(:,1:2) = true;
    else
        if useRidx
            if iscell(opts.Ridx)
                Ridxk = opts.Ridx{kk};
                ixfilt(Ridxk,1) = true;
            elseif isnumeric(opts.Ridx)
                Ridxk = opts.Ridx(kk,:);
                if isnan(Ridxk(1)); Ridxk(1) = 1; end
                if isnan(Ridxk(2)); Ridxk(2) = length(R); end
                ixfilt(Ridxk(1):Ridxk(2),1) = true;
            end
        end
        
        if opts.useRzLim
            RzIdxk = opts.RzLim(kk,:);
            if isnan(RzIdxk(1)); RzIdxk(1) = min(zR); end
            if isnan(RzIdxk(2)); RzIdxk(2) = max(zR); end
            ixfilt(:,2) = and( zR>=min(RzIdxk) , zR<=max(RzIdxk) );
        else
            RzIdxk = [];
        end  
        
        ixfilt = any(ixfilt,2) & sffilt;
        RidxAll = find( (any(ixfilt,2) & sffilt) );
    end
    
    % ------- (2) get dr/dz and z0 -------- 
    rStats = getRparams(zR,R,RidxAll); %,Ridxk);
    
    % Set z0 from fitlm or manual input
     if ~isempty(opts.z0)
         if size(opts.z0,1)==Ntracks
             z0k = opts.z0(kk,:);
         else
             z0k = opts.z0;
         end
         
         if size(z0k,2)==3
             z0  = z0k(2);
             z0ci = opts.z0(kk,[1 3]);
         elseif size(z0k,2)==1
             z0 = z0k;
             z0ci = [NaN NaN];
         end        
         if isnan(z0)
             useManualZ0 = false;
         else
             useManualZ0 = true;
         end
     else
         useManualZ0 = false;
     end
         
    if ~useManualZ0
        z0      = rStats.z0;
        z0ci    = rStats.z0ci;
    end
    
    % Get r0
    if isempty(opts.r0)
        r0      = rStats.r0;
        r0ci    = rStats.r0ci;
    elseif size(opts.r0,1)==Ntracks
        if size(opts.r0,2)==3
            r0   = opts.r0(kk,2);
            r0ci = opts.r0(kk,[1 3]);
        elseif size(opts.r0,2)==1
            r0   = opts.r0(kk);
            r0ci = [nan nan];
        end
    end
        
    
    % ------- (3) get T(z) -------- 
    
    [zT,T,sf,ff] = getTcurve(tk(kk),opts.Tstruct,opts.Tcurve,opts.useZmax);
    
    
   % Parse T indices/z limits and apply filter thresholds
    ixfilt = false(length(zT),2);
    sffilt = and(sf<=opts.satFracThresh, ff<=opts.filtFracThresh);
    if ~useTidx && ~opts.useTzLim
        ixfilt(:,1:2) = true;
    else
        if useTidx
            if iscell(opts.Tidx)
                Tidxk = opts.Tidx{kk};
                ixfilt(Tidxk,1) = true;
           elseif isnumeric(opts.Tidx)
                Tidxk = opts.Tidx(kk,:);
                if isnan(Tidxk(1)); Tidxk(1) = 1; end
                if isnan(Tidxk(2)); Tidxk(2) = length(T); end
%                 Tidxk = Tidxk(1):Tidxk(2);
                ixfilt(Tidxk(1):Tidxk(2),1) = true;
            end
        end
        
        if opts.useTzLim
            TzIdxk = opts.TzLim(kk,:);
            if isnan(TzIdxk(1)); TzIdxk(1) = min(zT); end
            if isnan(TzIdxk(2)); TzIdxk(2) = max(zT); end
            ixfilt(:,2) = and( zT>=min(TzIdxk) , zT<=max(TzIdxk) );
        end  
        
        ixfilt  = any(ixfilt,2) & sffilt;
        TidxAll = find(ixfilt);
    end
    
    
    
    % ------- (4) get Power law fit (B-value) -------- 
    if (zT(TidxAll(1))-z0)<0
        warning('Track %s-%i z0 estimate is unphyiscally high. Using z0=0',tk(kk).event,tk(kk).eventTrack)
        z0 = 0; % Could alternatively default to averaged image or a scaling value
        z0ci = [0 0]; 
        Rqc(kk) = 0;
    end
    switch opts.Bmethod
        case 'fitPL'
            Bfit = getTexponent(zT,T,TidxAll,...
                z0,z0ci,r0,opts.sensitivityFlag,opts.dT0Limits,opts.cLimits);
            if any(isnan(Bfit.Bpl))
                fprintf(' --> Missing B value: %s-%i\n',tk(kk).event,tk(kk).eventTrack)
                Tqc(kk)=0;
            end
    end

    
    % report/output opts (curves, idx , etc), zR,R,zT,T (full track versions),rStats,Bfit
    drdzci = rStats.Rmdl.coefCI(.01);
    
    trackFits(kk).eventName = tk(kk).event;
    trackFits(kk).trackNo   = tk(kk).eventTrack;
    trackFits(kk).Rqc       = Rqc(kk);
    trackFits(kk).Tqc       = Tqc(kk);
    trackFits(kk).zR        = zR;
    trackFits(kk).zT        = zT;
    trackFits(kk).R         = R;
    trackFits(kk).T         = T;
    trackFits(kk).Ridx      = RidxAll;
    trackFits(kk).Tidx      = TidxAll;
    trackFits(kk).rStats    = rStats;
    trackFits(kk).r0        = [r0ci(1) r0 r0ci(2)];
    trackFits(kk).z0        = [z0ci(1) z0 z0ci(2)];
    trackFits(kk).drdz      = [drdzci(2,1) rStats.Rmdl.Coefficients.Estimate(2) drdzci(2,2)];
    trackFits(kk).Bfit      = Bfit;
    trackFits(kk).Bpl       = Bfit.Bpl;

    fitArrays.tSpan(kk,:)   = tk(kk).t([1 end]);
    fitArrays.z0(kk,:)      = trackFits(kk).z0;
    fitArrays.r0(kk,:)      = trackFits(kk).r0;
    fitArrays.drdz(kk,:)    = trackFits(kk).drdz;
    fitArrays.r_r2(kk)      = trackFits(kk).rStats.Rmdl.Rsquared.Adjusted;
    fitArrays.r_rmse(kk)    = trackFits(kk).rStats.Rmdl.RMSE./mean(trackFits(kk).rStats.r); % Now normalized
    fitArrays.Bpl(kk,:)     = trackFits(kk).Bfit.Bpl;        % Power law fit (arguably only need this one)
    fitArrays.B_r2(kk)      = trackFits(kk).Bfit.gof{2}.rsquare;
    if ~and(0<=trackFits(kk).Bfit.gof{2}.rsquare,trackFits(kk).Bfit.gof{2}.rsquare<=1) % check for cases with failed fit
        fitArrays.B_nrmse(kk)   = nan;
    else
        fitArrays.B_nrmse(kk)   = trackFits(kk).Bfit.gof{2}.rmse;
    end
end



end