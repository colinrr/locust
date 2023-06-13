function [Vtrack,param] = trackVelocities(V,D,varargin) %param.iTrigger,param.trackWindow,param.detectWindow,param.Tpercentile)
% For any source window(s), with time indices given by win_tI and height
% indices given by src_zI, estimate vertical movement of the window using
% opticalFlow velocity vertical component. Track motion for all windows
% given by win_tI. 
%
%   V = velocity data cube struct from "thermOpticFlow.m"
%   D = thermal data cube struct from "getThermCube.m", with FIELDS:
%               T   : thermal data cube
%               M   : plume mask cube, dimensions equal to D.T 
%               z   : z position vector corresponding to 1st dimension of D.T
%               x   : x position vector corresponding to 2nd dimension of D.T
%               z   : z position vector corresponding to 1st dimension of D.T
%               t   : time position vector corresponding to 3rd dimension of D.T
%   param = pulse track parameter struct. See "pulseTrackInput.m"
%
%   param.iTriggerger   = starting frame indices for each track (ie detection frames)
%           --> scalar or vector of indices along 3rd DIM of Vz used
%           --> this is an output of thermPulseDetection.m
%   param.trackWindow = [I H] vertical position of TRACKING window.
%            --> I =  bottom pixel position
%            --> H =  height (in pixels) of tracking window.
%            Alternate: [Nx2] matrix, where N = length(param.iTriggerger), to use a
%            unique tracking window for each track
%
% OPTIONAL:
%   param.detectWindow = same as param.trackWindow, but for the DETECTION window if
%               thermPulseDetection was used. This provides an initial prior 
%               model for cluster detection.
%               --> Default = upper half of tracking window.
%               --> MIGHT CHANGE THIS TO ALLOW A FULL PRIOR SPECIFICATION
%   param.Tpercentile = Percentile value with which to filter out cool 
%                temperatures during clustering and tracking. 
%               --> Default = 70.
%
%
%
% OUT: V = struct containing velocities within tracking windows and
% calculated indices for statistics windows.
%
% I'll make this less confusing later...
%
% C Rowell Sep 2019
%  DEV NOTES:
%   1) Test with "verbose" output to refine objective f'n, then remove
%   unecessary outputs.
%      1a) "normal" vs "full" data output options for troubleshooting/testing
%   2) Implement backtracking from first window for complete pulse history
%   3) ? Use "findpeaks" on source spectra to estimate key timescales,
%   combine with 95th %ile detection velocities to get characteristic
%   lengthscales for a recommended tracking window size

%% VERSION HISTORY

% 1.0   Original tracker, designed and tested on event 25B, starting frames
%       137 - 140. 
%
% 1.01  Updates, tweaks, and revisions to get original demo tracks for 25B.
%       (See AGU2020 abstract). Original ensemble tests for best input 
%       values were run on this version. Last update ~July 30 2020.
%
% 1.10  A few key performance updates were made to allow more general
%       operation. Complete 25B and 24A tracks were achieved with this
%       version, finalized as of September 17 2020.
%
% 1.11  Initial runs for event 25A-4. Constrained maxClust for
%       initialization step (formerly was maxClust*2, now cannot be more
%       than input maxClust). Fixed a lower limit of 1 for pxTol due to
%       slow 25A velocities.
%           CHANGELOG:
%           --> Frame limiter fix to 'dzI' update during tracking window 
%               update.
%           --> Implemented a very light mask smoothing step in
%               restrictivePriorWarp to eliminate outlier pixels

clear textprogressbar
qcPlot = true;
%% PARSE
    % Get tracker
    param = loadPulseTrackInput(varargin{:});

    
    Nsources = length(param.iTrigger);
    
    disp('Grid calculations...')
    if isempty(param.trackWindow)
        nSrcWins = 1;
        nZ = round(sum(D.mask(1,:,param.iTrigger(1)))*0.5); % 1/2 width default
        param.trackWindow = [1 nZ];
    else
        nSrcWins = size(param.trackWindow,1); % Either 1 or same as length(param.iTrigger)
    end
    nZ = param.trackWindow(:,2); %size(src_zI,1);   % Height of initial statistics window(s)

    if isempty(param.detectWindow)
        param.detectWindow = param.trackWindow; % On the way to deprecating 2nd window - just set equal for now
%         param.detectWindow = [param.trackWindow(1)+round(param.trackWindow(2)/2) round(param.trackWindow(2)/2)];
    end
    if isempty(param.lambda)
        param.lambda = 0.1;
    end
    
    % Input checks
    % assert(and(isvector(src_zI),max(tI0)<length(D.t)),'src_zI needs to be a vector of indices in the range [1 size(win_tI,2)]')
    assert(and(all(param.iTrigger>=1), all(param.iTrigger<length(D.t))),...  
        '"param.iTrigger" indices must be between 1 and number of frames in "D"')
    assert(and(all(param.trackWindow(:)>=1), all(param.trackWindow(:)<= size(D.T,1))),...
        '"param.trackWindow" values must between 1 and vertical frame dimension')
    assert(or(nSrcWins==1,nSrcWins==Nsources),...
        'param.trackWindow must be either a [1x2] or [length(param.iTrigger)x2] matrix')
    assert(and(all(param.iTrigger>=1), all(param.iTrigger<length(D.t))),...
        '"param.iTrigger" indices must be between 1 and number of frames in "D"')

%     z = D.z;
%     dz = mean(diff(z));
    dt = mean(diff(V.t));
    Mf = length(D.z);
    Nf = length(D.x);
    
    % Calculate velocity tolerances
        
    if isempty(param.uMax) % A typical maximum velocity for any given frame
        param.uMax = prctile(max( sqrt(V.Vx.^2 + V.Vz.^2).*D.mask, [], [1 2]),96);
        fprintf('uMax:\t%.2f\n',param.uMax)
    end
    if isempty(param.uGrid) % Characteristic grid velocity
        param.uGrid = D.dx/mean(diff(D.t));
        fprintf('uGrid:\t%.2f\n',param.uGrid)
    end
    if isempty(param.pxTol) % Frame-to-frame motion tolerance (pixels)
        param.pxTol = max([param.uMax*param.uTol*dt/D.dz 1]); % Maximum allowable movement for a pixel
    end
    if isempty(param.memoryN)
        Vall = sqrt(V.Vz(D.mask).^2 + V.Vx(D.mask).^2);
        uVec = 0:0.1:max([param.uGrid prctile(Vall,99)]);
        V_pdf = zeros(length(uVec),10);
        for ii=1:size(V_pdf,2)
            Vsamp = randsample(Vall,1e4);
            V_pdf(:,ii) = ksdensity(Vsamp,uVec);
        end
        V_pdf = mean(V_pdf,2);
        [~,vI] = max(V_pdf);
        
        % 1/Modal Courant ratio
        param.memoryN = round(param.uGrid./uVec(vI));
        fprintf('memory:\t%i\n',param.memoryN)
    end
    
    % Window can always grow by at least one pixel, otherwise limited to uMax;
    growthSpeed = max([1 round(param.uMax/param.uGrid*param.uTol/2)]);
    
    % Stopping conditions
    manualStopTime = ~isempty(param.stopTime); 
    manualStopHeight = ~isempty(param.stopHeight);
    frameStop = false;
            
 %% Apply atmospheric profile removal?
    if isfield(D,'atmo') && param.applyAtmo
        % Get appropriate profile struct
        if iscell(D.atmo)
            if length(D.atmo)>1
                for aa=1:length(D.atmo)
                    Trms(aa) = D.atmo{kk}.Trms;
                end
                [~,atmoI] = min(Trms);
            else
                atmoI = 1;
            end
            atmo = D.atmo{atmoI};
        elseif isstruct(D.atmo)
            atmo = D.atmo;
        end
        
        zASL = D.z+D.z0ASL; % Image height vector ASL
        % Get atmospheric profile
        Tatmo = interp1(atmo.Height(~isnan(atmo.Height)),...
            atmo.Temperature(~isnan(atmo.Temperature)),zASL,'spline');

%         apply_atmo = true;
       
    else
        param.applyAtmo = false;
    end


if param.applyAtmo
    % Either remove atmo here or one frame at a time?
    threshT = -param.nHalfMax*(atmo.Tmode-atmo.T_halfMax);
    cax = [threshT 130]; %140];
else
    cax = [240 320];
end

%% Loop through each source time/frame
disp('Tracking features...')
for kk=1:Nsources
    Vtrack(kk).cubeDims = size(D.T);
    Vtrack(kk).N = size(V.Vz,3)-param.iTrigger(kk)+1; % Number of samples for this track
    fprintf('   Source %i/%i: %i frames\n',kk,Nsources, Vtrack(kk).N)
%     textprogressbar(sprintf('   Source %i/%i: ',kk,Nsources))
%     textprogressbar('')
    dispstat('','init')
    trackMaskTop = false;

    % SET UP: zI, detI(priors) for this track
    Vtrack(kk).Tprc = param.Tpercentile;
    Vtrack(kk).zI   = zeros([2 Vtrack(kk).N]);    
    Vtrack(kk).xI   = zeros([2 Vtrack(kk).N]);
    
     % Initialize vector of height indices
    if nSrcWins==1
        Vtrack(kk).zI(:,1) = [ param.trackWindow(1); sum(param.trackWindow)-1]; % Set source window
    elseif nSrcWins==Nsources
        Vtrack(kk).zI(:,1) = [param.trackWindow(kk,1); sum(param.trackWindow(kk,:))-1];
    end %--> OPTION HERE FOR FULLY MANUAL TRACK LOCATIONS
    
    % First: obtain velocities and integrate height for each frame
%     V(kk).Dzdt = zeros(V(kk).N,1); % Initialize local velocity
    Vtrack(kk).tI0  = param.iTrigger(kk);               % Time start index for this track
    Vtrack(kk).tI   = param.iTrigger(kk):numel(D.t);    % Time indices from original time vector
    Vtrack(kk).t    = V.t(param.iTrigger(kk):end);      % Get times for local velocity tracking
    
    % WINDOW
    Vtrack(kk).Tprc    = zeros(Vtrack(kk).N,1);         % Window percentiles
    Vtrack(kk).winIdx  = cell([Vtrack(kk).N 1]);        % Pixel indices for window
    Vtrack(kk).Vmu     = zeros(Vtrack(kk).N,1);         % Window mean W velocity
    Vtrack(kk).winZ    = zeros(Vtrack(kk).N,1);         % Window mean height vector
%     Vtrack(kk).winZ(1) = mean(z(Vtrack(kk).zI(:,1)));   % Initial centered window position
    
    % CLUSTER
%     Vtrack(kk).trackI  = zeros([nTZ Vtrack(kk).N]);   % Initialize vector of tracking indices
    Vtrack(kk).Nclust   = zeros(Vtrack(kk).N,1);        % Chosen number of clusters
%     Vtrack(kk).eigs     = zeros(15,Vtrack(kk).N);       % Spectral clustering eigenvalues
    Vtrack(kk).clustIdx = cell([Vtrack(kk).N 1]);       % Pixel indices for tracked cluster
    Vtrack(kk).clustI   = zeros(Vtrack(kk).N,3);         % Cluster mean, min, max I position
    Vtrack(kk).clustZ   = zeros(Vtrack(kk).N,1);         % Cluster mean Z position
    Vtrack(kk).clustT   = zeros(Vtrack(kk).N,1);         % mean tracked cluster Temperature
    Vtrack(kk).Tsum     = zeros(Vtrack(kk).N,1);         % mean tracked cluster Temperature    
    Vtrack(kk).clustW   = zeros(Vtrack(kk).N,1);         % mean tracked cluster W velocity
    Vtrack(kk).npx      = zeros(Vtrack(kk).N,1);         % number of pixels in window...?
    Vtrack(kk).iExtent  = zeros(Vtrack(kk).N,1);         % vertical extent of tracked cluster in px
    Vtrack(kk).jExtent  = zeros(Vtrack(kk).N,1);         % horizontal extent of tracked cluster in px
    Vtrack(kk).iCenter  = zeros(Vtrack(kk).N,1);         % i-coord for cluster tracking target
    Vtrack(kk).jCenter  = zeros(Vtrack(kk).N,1);         % j-coord for cluster tracking target
     
    
    % OBJECTIVE FUNCTION
%     Vtrack(kk).Cscore   = cell([Vtrack(kk).N 1]);       % Term 1 - data term
%     Vtrack(kk).Pvec     = cell([Vtrack(kk).N 1]);       % Prior vector - T, W, D, A scores
%     Vtrack(kk).Pscore   = cell([Vtrack(kk).N 1]);       % Term 2 - prior term (normed Pvec)
    Vtrack(kk).Obj      = cell([Vtrack(kk).N 1]);       % Final objective function: Cscore + lambda*Pscore

    % USEFUL TESTING STATS
    Vtrack(kk).Tprc      = zeros(Vtrack(kk).N,1); % Dynamic changes to Tprc
    Vtrack(kk).clustMedT = zeros(Vtrack(kk).N,1); % cluster median temp
    Vtrack(kk).clustVar  = zeros(Vtrack(kk).N,1); % cluster temp variance
    Vtrack(kk).winMedT   = zeros(Vtrack(kk).N,1); % window median temp
    Vtrack(kk).winVar    = zeros(Vtrack(kk).N,1);
%     Vtrack(kk).prcOver   = zeros(Vtrack(kk).N,1); % Percent overlap of intial cluster with prior

    
    trkPar       = param;
    % Set up prior memory
%     trkPar.priorMem.Nclust = zeros(trkPar.memoryN,1);
%     trkPar.priorMem.Tsum   = zeros(trkPar.memoryN,1);
%     trkPar.priorMem.Tbar   = zeros(trkPar.memoryN,1);
%     trkPar.priorMem.Wbar   = zeros(trkPar.memoryN,1);
%     trkPar.priorMem.Area   = zeros(trkPar.memoryN,1);
    trkPar.priorMem      = zeros(Mf,Nf,param.memoryN);
    trkPar.priorWarp     = cell(1,param.memoryN);
    trkPar.priorCenter   = zeros(Vtrack(kk).N,2); %zeros(param.memoryN,2);
    trkPar.priorExtent   = zeros(Vtrack(kk).N,2); %zeros(param.memoryN,2);
    trkPar.priorAspect   = zeros(Vtrack(kk).N,2); % zeros(param.memoryN,2);
    trkPar.priorMaxI     = zeros(Vtrack(kk).N,2);
    %     trkPar.priorMem.pvec   = zeros(param.memoryN,3); % Tsum,Wbar,Area
    
%     textprogressbar(' ')
    %% Loop through frames to get complete velocity set
    for jj=1:Vtrack(kk).N  % tI0(kk):size(Vz,3)
        dispstat(sprintf('    Frame %i/%i - %.1f%%',jj,Vtrack(kk).N,jj/Vtrack(kk).N*100))
        
        % Get data for this frame
        Tj    = D.T(:,:,(jj-1)+param.iTrigger(kk));
        Wj    = V.Vz(:,:,(jj-1)+param.iTrigger(kk));
        Uj    = V.Vx(:,:,(jj-1)+param.iTrigger(kk));
%         zIj   = Vtrack(kk).zI(:,jj);              % Get current stats window
        maskj = D.mask(:,:,(jj-1)+param.iTrigger(kk));      % Plume mask
%         dTdzj = dTdz(:,:,(jj-1)+param.iTrigger(kk));   % Get vertical temperature gradient   

        % Atmo profile? Better here or before loop?
        if param.applyAtmo
            %    T - profile - DC-shift
            Tj = Tj-Tatmo-atmo.T_halfMax;
            threshI = Tj<threshT;
            Tj(threshI) = threshT;
            maskj(threshI) = 0;
        end
        
        %% Tracking window: retrieve velocities
        
        % Params for this time step
%        trkPar       = param;
       trkPar.zI    = Vtrack(kk).zI(:,jj);
       trkPar.xI    = Vtrack(kk).xI(:,jj);
%         trkPar.prior = prior;
        

%% INITIALIZATION STEP ON FIRST FRAME
          % --> WARP cluster mask based on velocity fields, round nearest to get new
        % mask
        
        % Get/reset initial window positioning
        if jj==1
            

            
            % Check for plume mask top
            if param.maskTracking
                masktop = find(any(maskj,2),1,'last');
                if and(masktop>=trkPar.zI(1),masktop<=trkPar.zI(2))
                    trackMaskTop = true;
%                     trkPar.lambda = 0.05;
                end
            end
            
            if trackMaskTop
                trkPar.zI(2) = masktop;
                trkPar.xI = [find(any(maskj(trkPar.zI(1):trkPar.zI(2),:),1),1,'first');...
                             find(any(maskj(trkPar.zI(1):trkPar.zI(2),:),1),1,'last')];
%                 dzI = [pixBuffer-diff(trkPar.zI); pixBuffer];
%                 toproi = getROI(D.mask(:,:,9),'iLims',trkPar.zI);
%                 trkPar.xI = toproi(3:4)';
                
                Vtrack(kk).zI(:,jj) = trkPar.zI;
                Vtrack(kk).xI(:,jj) = trkPar.xI;
                winAspect = (range(trkPar.zI)+1)./(range(trkPar.xI)+1);
                
%                 iExtent0 = range(toproi(1:2));
%                 jExtent0 = range(toproi(3:4));

%                 if isempty(param.winSzRatio)
% %                     trkPar.winSzRatio = 1.5; %param.winSzRatio;
%                     trkPar.winSzRatio = param.trackWindow(2)/iExtent0; %param.winSzRatio;
%                 end

                % In mask tracking mode, we want large clusters close to
                % mask scale
                switch trkPar.clusterMode
                    case 'full'
                        trkPar.minClust = 1;
                        trkPar.maxCust  = 2;
                    case 'fast'
                        trkPar.minClust = 2;
                        trkPar.maxCust  = 4;                        
                    case 'eig'
                        trkPar.minClust = 2;
                        trkPar.maxCust  = 8;
                end
                
            else
            
                % Run initial
                trkPar0 = trkPar;
                trkPar0.clusterMode = 'eig';
                trkPar0.maxClust = min([trkPar0.maxClust 12]); % 12 is a realistic max
                S = pickVelocities(Tj,maskj,Wj,Uj,D.x,D.z,trkPar0);

                [subI,subJ]  = ind2sub(size(maskj),S.trackedIdx);
                iExtent0 = range(subI);
                jExtent0 = range(subJ);

                % Set window relative positioning - minimum buffer of 3 pixels
                pixBuffer = max([3 ceil(param.uTol.*(S.regionAvgs(S.trackedCluster,5)/D.dx*dt))]);
                dzI = [pixBuffer-diff(trkPar.zI); pixBuffer];
    %             dzI = trkPar.zI-max(subI); 

                % Cutoff positions for cluster leading edge to trigger
                % out-of-frame tracking stop: Keep detection window in frame
                minIwarp = max(subI)-trkPar.detectWindow(1)+1; 
                maxIwarp = Mf;

                % Re-do initial window scaling... and run again?
%                 if isempty(param.winSzRatio)
    %                 param.winSzRatio = param.trackWindow(2)/iExtent0;
%                     trkPar.winSzRatio = 1.5; %param.winSzRatio;
    %                 trkPar.winSzRatio = param.trackWindow(2)/iExtent0; %param.winSzRatio;
%                 end
                winAspect = min([1 iExtent0/jExtent0]); % Window restriction - at least as wide as tall
                Vtrack(kk).jCenter(jj)   = round(mean(subJ));

                % update window sizing/positioning to run again
                dzI       = [dzI(2)-round(iExtent0.*trkPar.winSzRatio)+1; dzI(2)];
                dxI       = [-1; 1].*ceil((diff(dzI)+1)./winAspect*0.5);

                % Set up to run again
                % -> adjust Tprc for similar amount of input data...
                % -> adjust window size with a regression
    %             Vtrack(kk).iCenter(jj) = Vtrack(kk).zI(2,jj) - dzI(2);
                Vtrack(kk).iCenter(jj) = max(subI);
                Vtrack(kk).zI(:,jj) = Vtrack(kk).iCenter(jj)+dzI;
                Vtrack(kk).xI(:,jj) = Vtrack(kk).jCenter(jj)+dxI;
    %             S = pickVelocities(Tj,maskj,Wj,Uj,D.x,D.z,trkPar);

                % Check for frame edges
                Vtrack(kk).zI(Vtrack(kk).zI(:,jj)<1,jj) = 1;
                Vtrack(kk).zI(Vtrack(kk).zI(:,jj)>Mf,jj) = Mf;
                Vtrack(kk).xI(Vtrack(kk).xI(:,jj)<1,jj) = 1;
                Vtrack(kk).xI(Vtrack(kk).xI(:,jj)>Nf,jj) = Nf;

                trkPar.zI    = Vtrack(kk).zI(:,jj);
                trkPar.xI    = Vtrack(kk).xI(:,jj);
                
                % Re-calculate number of clusters - pretty touchy at
                % the moment...
                Awin = (trkPar.winSzRatio*iExtent0).^2/winAspect; % Optimimum window size
                NclustOpt = Awin/numel(S.trackedIdx); % roughly optimum number of clusters
                
                % Alt formula
                Awin = (diff(trkPar.xI)+1)*(diff(trkPar.zI)+1); % Approximate window area
                NclustOpt = (Awin*(1-trkPar.Tpercentile/100))/numel(S.trackedIdx); % roughly optimum number of clusters

                switch trkPar.clusterMode
                    case 'full'
%                         trkPar.minClust  = min([max([trkPar.minClust round(NclustOpt)]) trkPar.maxClust-1]);
%                         trkPar.maxClust  = min([min([trkPar.minClust+3 ceil(1.75*NclustOpt)]) trkPar.maxClust]);
                        
                        trkPar.minClust = min([max([param.minClust round(NclustOpt)-1]) param.maxClust-1]);
                        trkPar.maxClust  = max([min([param.maxClust round(NclustOpt)+2]) param.minClust+1]);
                        
%                     case 'fast'
%                         
%                         trkPar.minClust = max([trkPar.minClust round(NclustOpt)]);
%                         trkPar.maxClust  = min([trkPar.minClust+2 ceil(1.75*NclustOpt)]);
%                         error
%                     case 'eig'
%                         trkPar.minClust = max([trkPar.minClust round(NclustOpt)]);
%                         trkPar.maxClust  = min([trkPar.minClust+6 ceil(1.75*NclustOpt)]);
%                         error
                    otherwise
                        error
    
                end

            end
        end     
        
%% Do the thing

        S = pickVelocities(Tj,maskj,Wj,Uj,D.x,D.z,trkPar);
        cI = S.trackedCluster;
        % Current cluster
        [subIjj,subJjj]  = ind2sub(size(maskj),S.trackedIdx);
        
%% RECORD TRACKING DATA...
          %--> trackedIdx -> cluster
          %--> tracking window - new position based on warped prior?
        % WINDOW
        Vtrack(kk).Tprc(jj)    = S.Tpercentile;                  % Varying percentile?
        Vtrack(kk).winIdx{jj}  = S.winIdx;                       % Pixel indices for window
        Vtrack(kk).Vmu(jj)     = mean(Wj(S.winIdx));             % Window mean W velocity
        Vtrack(kk).winZ(jj)    = mean(D.z(Vtrack(kk).zI(:,jj))); % Window mean height vector
        
        % CLUSTER
        Vtrack(kk).Nclust(jj)       = S.NclustCalc;
%         Vtrack(kk).eigs(:,jj)   = round(S.eigd,8);
        Vtrack(kk).clustIdx{jj} = S.trackedIdx;               % Cluster indices
%         Vtrack(kk).clustZ(jj)   = S.regionAvgs(cI,2);         % Cluster mean Z position
        Vtrack(kk).clustT(jj)   = S.regionAvgs(cI,3);         % mean tracked cluster Temperature
        Vtrack(kk).Tsum(jj)     = sum(Tj(S.trackedIdx));      % sum tracked cluster Temperature  
        Vtrack(kk).clustW(jj)   = S.regionAvgs(cI,5);         % mean tracked cluster W velocity 
        Vtrack(kk).clustI(jj,:) = [min(subIjj) mean(subIjj) max(subIjj)];  % Cluster I positioning
        Vtrack(kk).npx(jj)      = S.clusterLargestRegion(cI); % number of pixels in cluster

        % OBJECTIVE FUNCTION
        Vtrack(kk).Obj{jj}      = S.Obj; % = [Nclust x 8] ; cols = [clusterIndex Cscore P_T P_W P_D P_A Pscore Objective_score]

        % USEFUL TESTING STATS
        Vtrack(kk).clustMedT(jj) = median(Tj(S.trackedIdx)); % cluster median temp
        Vtrack(kk).clustVar(jj)  = var(Tj(S.trackedIdx)); % cluster temp variance
        Vtrack(kk).winMedT(jj)   = median(Tj(S.winIdx)); % window median temp
        Vtrack(kk).winVar(jj)    = var(Tj(S.winIdx));
%         if jj>1
%             Vtrack(kk).prcOver(jj)   = numel(intersect(S.rawTrackedIdx,find(trkPar.prior.mask)))./numel(S.rawTrackedIdx)*100; % Percent overlap of intial cluster with prior
%         end
%% UPDATE PRIOR

        % UPDATE PRIOR MASK MEMORY        
        % ---- (A) Warp CURRENT cluster for prior -------
%         subIwarp     = round(subIjj + Wj(S.trackedIdx)./V.dx*dt);
%         subJwarp     = round(subJjj + Uj(S.trackedIdx)./V.dz*dt);
%         Idxwarp      = sub2ind(size(maskj),subIwarp,subJwarp);
%         warpMask     = false(size(maskj));
%         warpMask(Idxwarp) = true;
% %         warpMask     = imfill(warpMask);
%         trkPar.priorMem = cat(3, trkPar.priorMem(:,:,2:end), warpMask);
%         priorIn = squeeze(any(trkPar.priorMem,[1 2]));
        % ----- (B) Warp ALL priors ------
        trkPar.priorWarp = [trkPar.priorWarp(2:end) {[subIjj subJjj]}];
        priorIn = cell2mat(cellfun(@(x) any(x(:)),trkPar.priorWarp,'UniformOutput',false));
%         trkPar.priorMem = cat(3, trkPar.priorMem(:,:,2:end), warpMask);
        for pi = find(priorIn)
%             [subIpi,subJpi]  = ind2sub(size(maskj),trkPar.priorMem(:,:,pm));
            Idxwarp      = sub2ind(size(maskj),round(trkPar.priorWarp{pi}(:,1)),round(trkPar.priorWarp{pi}(:,2))); % First index
            trkPar.priorWarp{pi}(:,1)     = trkPar.priorWarp{pi}(:,1) + Wj(Idxwarp)./V.dx*dt;
            trkPar.priorWarp{pi}(:,2)     = trkPar.priorWarp{pi}(:,2) + Uj(Idxwarp)./V.dz*dt;
            
            % Clip any pixels that warp out of frame
            iIn = and(round(trkPar.priorWarp{pi}(:,1))>=1,round(trkPar.priorWarp{pi}(:,1))<=Mf);
            ijIn = iIn & and(round(trkPar.priorWarp{pi}(:,2))>=1,round(trkPar.priorWarp{pi}(:,2))<=Nf);
            if any(~ijIn)
                trkPar.priorWarp{pi} = trkPar.priorWarp{pi}(ijIn,:);
            end
            
            Idxwarp      = sub2ind(size(maskj),round(trkPar.priorWarp{pi}(:,1)),round(trkPar.priorWarp{pi}(:,2))); % Warped index
            warpMask     = false(size(maskj));
            warpMask(Idxwarp) = true;
            
            % Fill any holes
%             warpImg = regionprops(warpMask,'FilledImage','BoundingBox','Area');
%             bbox = [ceil(warpImg.BoundingBox(1:2)) floor(sum(warpImg.BoundingBox([1 3]))) floor(sum(warpImg.BoundingBox([2 4])))];
%             warpMask(bbox(2):bbox(4),bbox(1):bbox(3)) = warpImg.FilledImage;
            
            % Update prior
            trkPar.priorMem(:,:,pi) = warpMask; % Update prior mask
        end
        subIwarp = trkPar.priorWarp{end}(:,1);
        subJwarp = trkPar.priorWarp{end}(:,2);
        % ------------------------
        

        if jj==1
%             Idxwarp = S.trackedIdx;
            [subI,subJ]  = ind2sub(size(maskj),S.trackedIdx);
        else
%             Idxwarp = find(sum(trkPar.priorMem,3)>=2);
            [subI,subJ]  = ind2sub(size(maskj),find(sum(trkPar.priorMem,3)>=2));
        end
%         trkPar.priorCenter = cat(1, trkPar.priorCenter(2:end,:), [max(subIwarp) round(mean(subJwarp))]);
%         trkPar.priorExtent = cat(1, trkPar.priorExtent(2:end,:), [range(subIwarp) range(subJwarp)]);
%         trkPar.priorCenter = cat(1, trkPar.priorCenter(2:end,:), [max(subI) round(mean(subJ))]);
%         trkPar.priorExtent = cat(1, trkPar.priorExtent(2:end,:), [range(subI) range(subJ)]);
        trkPar.priorCenter(jj,:) = [max(subI) mean(subJ)];
        trkPar.priorExtent(jj,:) = [range(subI) range(subJ)];
        trkPar.priorMaxI(jj)     = max(subIjj);
        
        
        % BUILD NEXT PRIOR...
    % Mean or weighted mean?
        jjIn = (jj-param.memoryN+1):jj;
        priorStartI = max([jjIn(1) 1]);
        trkPar.prior.Nclust = S.NclustCalc;
        trkPar.prior.Tsum   = mean(Vtrack(kk).Tsum(priorStartI:jj));
        trkPar.prior.Wbar   = mean(Vtrack(kk).clustW(priorStartI:jj));
        trkPar.prior.Area   = mean(Vtrack(kk).npx(priorStartI:jj));
%         trkPar.prior.Trange = [prctile(Tj(S.trackedIdx),2.5) prctile(Tj(S.trackedIdx),98.5)];
%         trkPar.prior.Trange = [min(Tj(S.trackedIdx)) max(Tj(S.trackedIdx))];
        trkPar.prior.Trange = [prctile(Tj(S.trackedIdx),2.5) max(Tj(S.trackedIdx))];
        
        keepI = min([jj,3,trkPar.memoryN]); % Min # of prior timesteps for keeping pixels
        
%         if trackMaskTop % Experimental: Force prior to track with plume front during trackMaskTop mode
%             [~,trkPar.prior.mask] = getROI(maskj,'iLims',Vtrack.zI(:,jj)+[round(diff(Vtrack.zI(:,jj))/2); 0]);
%         else
            
          trkPar.prior.mask   = double(sum(trkPar.priorMem,3)>=keepI);
            
%         end
%         if jj==1
%             trkPar.prior.mask   = double(sum(trkPar.priorMem,3)>0); % Keep first mask
%         else
%             trkPar.prior.mask   = double(sum(trkPar.priorMem,3)>1); % Keep pixels showing up at least twice
%         end

        
        % Extents from tracked cluster
%         Vtrack(kk).iExtent(jj) = range(subI);
%         Vtrack(kk).jExtent(jj) = range(subJ);
        % Extents from most recent prior
%         meanPrior = round(mean(trkPar.priorMem(:,:,any(trkPar.priorMem,[1 2])),3));
%         Vtrack(kk).iExtent(jj) = sum(any(trkPar.priorMem(:,:,end),2),1);
%         Vtrack(kk).jExtent(jj) = sum(any(trkPar.priorMem(:,:,end),1),2);
        % Averaged prior extents
%         priorIn = squeeze(any(trkPar.priorMem,[1 2]));
        Vtrack(kk).iExtent(jj) = mean(sum(any(trkPar.priorMem(:,:,priorIn),2),1));
        Vtrack(kk).jExtent(jj) = mean(sum(any(trkPar.priorMem(:,:,priorIn),1),2));
        % --> Characteristic aspect ratio
        trkPar.priorAspect(jj) = min([ 1 mean(Vtrack(kk).iExtent(1:jj)./Vtrack(kk).jExtent(1:jj))]);

        %% UPDATE TRACKING WINDOW        

        % Get cluster target pixel (cluster top-center for now)
%         Vtrack(kk).iCenter(jj)   = Vtrack(kk).zI(2,jj) - dzI(2);
        Vtrack(kk).iCenter(jj)   = mean(trkPar.priorMaxI(jjIn(priorIn)));
        Vtrack(kk).jCenter(jj)   = mean(subJ);
        Vtrack(kk).clustZ(jj)   = (Vtrack(kk).iCenter(jj)-1)*D.dz+D.z(1); % Cluster tracked Z position (approx leading edge)

        if trackMaskTop
%             topNextMask = find(any(D.mask(:,:,(jj)+param.iTrigger(kk)),2),1,'last');
            nextROI = getROI(D.mask(:,:,(jj)+param.iTrigger(kk)),'iLims',[trkPar.zI(1) length(D.z)]);
%             Vtrack(kk).xI(:,jj+1) = nextROI(3:4)';
            
            % Original xI lims setting for track window
%             Jsz = diff(Vtrack(kk).xI(:,jj))+1;
%             targetJsz = range(nextROI(3:4))+1;
%             diffJsign = sign(targetJsz - (diff(Vtrack(kk).xI(:,jj))+1));
%             diffJsz = diffJsign.*round(min(abs([growthSpeed targetJsz - (diff(Vtrack(kk).xI(:,jj))+1)])));  
%             Vtrack(kk).xI(:,jj+1) = round(mean(Vtrack(kk).xI(:,jj)) + ((Jsz+diffJsz)).*[-1/2 1/2]);
             
            % Alt xI lims setting for track window
            diffJ = nextROI(3:4)' - Vtrack(kk).xI(:,jj);
            diffJsign = sign(diffJ);
            diffJ = min(abs([diffJ growthSpeed.*[1;1]]),[],2).*diffJsign;
            Vtrack(kk).xI(:,jj+1) = Vtrack(kk).xI(:,jj) + diffJ;
            
            targetIsz = (diff(Vtrack(kk).xI(:,jj+1))+1).*winAspect;
            
            Isz = diff(Vtrack(kk).zI(:,jj))+1;
            diffIsign = sign(targetIsz - (diff(Vtrack(kk).zI(:,jj))+1));
            diffIsz = diffIsign.*round(min(abs([growthSpeed targetIsz - (diff(Vtrack(kk).zI(:,jj))+1)])));
            
            Vtrack(kk).zI(:,jj+1) = [nextROI(2)-(Isz + diffIsz)+1;...
                                     nextROI(2)];
            dzI = round(trkPar.zI - Vtrack.iCenter(jj));
            dxI = round(trkPar.xI - Vtrack.jCenter(jj));
                                 
        else
            nextiCenter  = max(subIwarp);
            nextjCenter  = mean(subJwarp);

            % limit window tracking motion with pxTol
            diCenter = min(abs([nextiCenter - Vtrack(kk).iCenter(jj) floor(param.pxTol)]));
            diSign = sign(nextiCenter - Vtrack(kk).iCenter(jj));
    %         djCenter = min(abs([round(nextjCenter - Vtrack(kk).jCenter(jj)) floor(param.pxTol)]));
            djCenter = min(abs([nextjCenter - Vtrack(kk).jCenter(jj) floor(param.pxTol)]));
            djSign = sign(nextjCenter - Vtrack(kk).jCenter(jj));

            % Aspect ratio/dimensions
            % --> Characteristic cluster dimension
                % update window sizing/positioning to run again
    %         targetIsz       = round(mean(Vtrack(kk).iExtent(max([1 jj-param.memoryN+1]):jj).*trkPar.winSzRatio));
            targetIsz       = Vtrack(kk).iExtent(jj).*trkPar.winSzRatio;
            targetJsz       = max([Vtrack(kk).jExtent(jj)*1.25 ceil(targetIsz./trkPar.priorAspect(jj))]); % From current prior
    %         targetJsz       = ceil((diff(dzI)+1)./winAspect);       % From original aspect

            % Signs and sizes of window dimension change
            diffIsign = sign(targetIsz - (diff(Vtrack(kk).zI(:,jj))+1));
            diffJsign = sign(targetJsz - (diff(Vtrack(kk).xI(:,jj))+1));
            diffIsz = diffIsign.*round(min(abs([growthSpeed targetIsz - (diff(Vtrack(kk).zI(:,jj))+1)])));
            diffJsz = diffJsign.*round(min(abs([growthSpeed targetJsz - (diff(Vtrack(kk).xI(:,jj))+1)])));  

            Jsz = diff(Vtrack(kk).xI(:,jj))+1 + diffJsz;

            % J size adjustment directionailty

            % Currently keeps window top at fixed position relative to cluster
            % top/"front" - some vertical directionality is implied. Could be
            % generalized by using vector velocity of tracked cluster.
            dzI = [max([dzI(1) - diffIsz -(round(Vtrack(kk).iCenter(jj) + diSign.*diCenter)-1)]) dzI(2)];


            % Update tracking window X position - Stays centered in X/J
            if mod(Jsz,2)==0
                Usign   = sign(mean(Uj(S.trackedIdx)));
                if Usign>=0
                    dxI = Jsz/2*[-1 1] + [1 0];
                else
                    dxI = Jsz/2*[-1 1] + [0 -1];
                end
            else
                dxI = fix(Jsz/2*[-1 1]);
            end

    %         Vtrack(kk).zI(:,jj+1) = round(Vtrack(kk).iCenter(jj)) + dzI + diSign.*diCenter;
    %         Vtrack(kk).xI(:,jj+1) = round(Vtrack(kk).jCenter(jj)) + dxI + djSign.*djCenter;
            Vtrack(kk).zI(:,jj+1) = round(Vtrack(kk).iCenter(jj) + diSign.*diCenter) + dzI;
            Vtrack(kk).xI(:,jj+1) = round(Vtrack(kk).jCenter(jj) + djSign.*djCenter) + dxI;
        end
        % Check for frame edges
        Vtrack(kk).zI(Vtrack(kk).zI(:,jj+1)<1,jj+1) = 1;
        Vtrack(kk).zI(Vtrack(kk).zI(:,jj+1)>Mf,jj+1) = Mf;
        Vtrack(kk).xI(Vtrack(kk).xI(:,jj+1)<1,jj+1) = 1;
        Vtrack(kk).xI(Vtrack(kk).xI(:,jj+1)>Nf,jj+1) = Nf;
        
        % POSSIBLY BACKTRACK TO PREVIOUS FRAME TO GET FEATURE INITIAL EVOLUTION?
        
 %%       
%         % Get window indices centered around new z
%         [~,Vtrack(kk).zI(:,jj+1)] = getMinK(abs(z-Vtrack(kk).z(jj+1)), nZ);
        
%         textprogressbar(jj/Vtrack(kk).N*100)

%% Temporary plotting for QC'ing

    if qcPlot
        if jj==1
            tempFig = figure('position',[200 100 1000 1000]);
        end
        cmask = false(size(maskj));
        rcmask = cmask;
        cmask(Vtrack(kk).clustIdx{jj}) = true;
        rcmask(S.rawTrackedIdx) = true;
        [~,~,cpoly] = getROI(cmask,'maxRegions',1);
        [~,~,rcpoly] = getROI(rcmask,'maxRegions',1);        
        [~,~,ppoly] = getROI(sum(trkPar.priorMem,3)>1,'maxRegions',1);
        [~,~,wpoly] = getROI(maskj,'iLims',Vtrack(kk).zI(:,jj),'jLims',Vtrack(kk).xI(:,jj),'maxRegions',1);
        [izoom,~,~] = getROI(maskj,'pad',20);
%         figure(tempFig)
        clf(tempFig)
        pax = gca;
        imagesc(pax,Tj)
        colormap(pax,gray(250))
        caxis(pax,cax)
        colorbar(pax)
        set(pax,'YDir','normal');
        hold on
        pl(1) = plot(pax,wpoly.X,wpoly.Y,'g','LineWidth',2);
        pl(2) = plot(pax,rcpoly.X,rcpoly.Y,'m','LineWidth',1);
        if ~isempty(ppoly)
            pl(3) = plot(pax,ppoly.X,ppoly.Y,'c','LineWidth',1);
            pl(4) = plot(pax,cpoly.X,cpoly.Y,'b','LineWidth',2);
            llab  = {'Tracking Window','Initial Cluster','Prior Memory','Warped Cluster'};
        else
            llab  = {'Tracking Window','Initial Cluster','Warped Cluster'};
            pl(3) = plot(pax,cpoly.X,cpoly.Y,'b','LineWidth',2);
        end
%         axis([20 360 0 450])
        title(pax,sprintf('jj = %i, Idx = %i, t = %.2f s',jj,(jj-1)+param.iTrigger(kk),V.t((jj-1)+param.iTrigger(kk))))
        legend(pax,pl,llab)
        %         xlim([izoom(3:4)])
%         ylim([izoom(1:2)])
%         axis([40 250 1 300])
        axis(pax,[1 300 1 250])
        pause(0.05)
    end  
        
%         [~,levelMask] = getROI(maskj,'iLims',Vtrack.zI(:,jj));
%         Tlevel = Tj(levelMask);
%         edges = min(Tlevel):1:max(Tlevel);
%         medLevel(jj) = median(Tlevel);
%         stdLevel(jj) = std(Tlevel);
%         medClust(jj) = median(Tj(S.trackedIdx));
%         stdClust(jj) = std(Tj(S.trackedIdx));
%         
        % QC stop
        if jj>=70
            disp('Flargh!')
%             fprintf('%s\n',mat2str(dzI))
        end
        
%% Stopping conditions

        % Condition for tracking into the frame bottom/top
%         if or(nextiCenter>maxIwarp,nextiCenter<minIwarp)
% %             Vtrack(kk).N = jj-1;
%             % Truncate structs here
%             Vtrack(kk) = truncateVtrack(Vtrack(kk),jj);
%             disp('Track ended: tracked beyond image bounds.')
%             break
% %         elseif Iwarp<minIwarp
%         end        

        % 1: Out of frames
        
        % TEMPRORARY STOP CONDITION FOR TESTING
        if manualStopTime
            if Vtrack.t(jj)>=param.stopTime
                frameStop = true;
            end
        end
        if manualStopHeight
            if Vtrack.clustZ(jj)>=param.stopHeight
                frameStop = true;
            end
        end
%         if D.z(round(Vtrack(kk).iCenter(jj)))>700   % Height condition
        % Vtrack.t(jj)>=58 %90 % 20 % 55               % Time condition
        if or(frameStop,jj==Vtrack(kk).N) 
            Vtrack(kk) = truncateVtrack(Vtrack(kk),jj);
            Vtrack(kk).stopCode = 1;
            disp('Track ended: Out of frames.')
            break
        end

        % 2: Tracking off frame edge: window is at frame edge and has size close
        % to cluster size
        if or(and(any(or(Vtrack(kk).zI==1,Vtrack(kk).zI==Mf)),abs(diff(dzI)+1 - Vtrack(kk).iExtent(jj))/Vtrack(kk).iExtent(jj)<0.1),...
              and(any(or(Vtrack(kk).xI==1,Vtrack(kk).zI==Nf)),abs(diff(dxI)+1 - Vtrack(kk).jExtent(jj))/Vtrack(kk).jExtent(jj)<0.1))
            % Truncate structs here
            Vtrack(kk) = truncateVtrack(Vtrack(kk),jj);
            Vtrack(kk).stopCode = 2;
            disp('Track ended: tracked beyond image bounds.')
            break
        end
        
        % 3: Cluster is well-mixed, statistically indistinguishable
        
        % 4: Tracking has wandered off to a different feature
        if S.offTrack
            Vtrack(kk) = truncateVtrack(Vtrack(kk),jj);
            Vtrack(kk).stopCode = 4;
            disp('Track ended: detected cluster is off track.')
            break
        end        
    end % end loop over frames
    Vtrack(kk).ObjLabels = {'clusterIdx','Cscore','T','W','D','A','Pscore','Obj'};

%     disp('check in...')
%     textprogressbar(' -> Done')

    %% Objective function plots
%     numClusters = 2;
%     trackLength = jj; %Vtrack(kk).N
%     objScore = zeros(numClusters,size(S.Obj,2),trackLength);
%     for ll=1:trackLength
%         objScore(:,:,ll) = Vtrack(kk).Obj{ll}(1:numClusters,:);
%     end
%     nClustVec = Vtrack(kk).Nclust(1:trackLength);
%     objScore = permute(objScore,[3 1 2]);
%     
%     dyp=0.02;
%     figure('name',sprintf('Track %i',kk),'position',[100 100 800 1200])
%     tightSubplot(8,1,1,[],dyp)
%     plot(nClustVec) % Nclust
%     ylabel('# Clusters')
%     set(gca,'FontSize',12)
%     title(sprintf('Objective function terms for 3 best clusters, lambda = %.1f',param.lambda))
%     tightSubplot(8,1,2,[],dyp)
%     plot(objScore(:,:,2)) % Cscore
%     ylabel('Cscore')
%     set(gca,'FontSize',12)
%     tightSubplot(8,1,3,[],dyp)
%     plot(objScore(:,:,7)) % Pscore
%     ylabel('Pscore')
%     set(gca,'FontSize',12)
%     tightSubplot(8,1,4,[],dyp)
%     plot(objScore(:,:,8),'LineWidth',2) % ObjFun
%     ylabel('OBJ FUN')
%     set(gca,'FontSize',12)
%     tightSubplot(8,1,5,[],dyp)
%     plot(objScore(:,:,3)) % T
%     ylabel('P_T')
%     set(gca,'FontSize',12)
%     tightSubplot(8,1,6,[],dyp)
%     plot(objScore(:,:,4)) % W
%     ylabel('P_W')
%     set(gca,'FontSize',12)
%     tightSubplot(8,1,7,[],dyp)
%     plot(objScore(:,:,5)) % D
%     ylabel('P_D')
%     set(gca,'FontSize',12)
%     tightSubplot(8,1,8,[],dyp)
%     plot(objScore(:,:,6)) % A
%     ylabel('P_A')
%     xlabel('Frame #')
%     set(gca,'FontSize',12)

    %% old crap
    % Smooth out the mess?
%     V(kk).Vsmooth = smooth(V(kk).Vmu,20);
    % Integrate to get height
%     V(kk).z = cumtrapz(V(kk).t,V(kk).Vsmooth);
    
    % After acquiring velocity values, will smooth, integrate and discretize to
    % obtain window height indices
%     V(kk).zI = zeros(size(src_zI,1),V(kk).N);
%     V(kk).zI(:,1) = src_zI;
%     
%     
%     [zdiff,midx] = sort(abs(V(kk).z - z));
    

%     V(kk).NumWins = size(src_tI(:,tI0(kk):end),2);
    
%     V(kk).winZ = zeros(V(kk).NumWins,1);
%     V(kk).winZ(1) = mean(z(src_zI));
%     V(kk).wint = mean(t(src_tI(:,tI0(kk):end)),1)'; % Window centered times
    
%     Vup = 0;
%     Vcut = zeros(NumVals,S.NumWins);


%     zI = src_zI;
%     for jj=tI0(kk):V(kk).N-1
%         
%         % Get values from minicube...or only one frame?
%         Vcut = Vz(zI,:,tI0(kk)) .* M(zI,:,src_tI(:,jj));
%         Vcut      = Vcut(Vcut~=0);
% 
%         % Get "mean" upwards flow velocity
%         Vup = mean(Vcut);  % Couple different measures I should try here
%         
%         Vup_px = round(Vup/dz*dt);
%         % Find new zI
%         zI = zI+Vup_px;
%        
%                 % Record current window, get masked Vcube
%         V(kk).zI(:,jj+1) = zI;
%         V(kk).winZ(jj+1) = V(kk).winZ+Vup_px*dz;
%     end
end


end

% function Vtrack  = truncateVtrack(Vtrack,N)
% 
%     ff = fieldnames(Vtrack);
%     for fi = 1:length(ff)
%         fn = ff{fi};
%         fSz = size(Vtrack.(fn));
%         dimChk = fSz==Vtrack.N;
%         if any(dimChk)
%             for dd=1:ndims(Vtrack.(fn))
%                 if dimChk(dd)
%                     idx{dd} = 1:N;
%                 else
%                     idx{dd} = 1:size(Vtrack.(fn),dd);
%                 end
%             end
% %             if iscell(Vtrack.(fn))
% %                 Vtrack.(fn) = Vtrack.(fn){idx{:}};
% %             elseif isnumeric(Vtrack.(fn))
%             Vtrack.(fn) = Vtrack.(fn)(idx{:});
% %             end
%         end
%     end
%     Vtrack.N = N;
% end

function [vals,idx] = getMinK(A,n)
% get n smallest elements in array A
[Asorted, Aidx] = sort(A(:));
vals = Asorted(1:n);
idx  = sort(Aidx(1:n));
end