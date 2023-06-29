function [S] = getCluster(T,mask,W,U,x,z,param)
% [Vo,To,pixIdx,Vmu,roi] = getCluster(T,W,mask,U,x,z,varargin)
% REQUIRED:
%       T           = single thermal image frame
%       mask        = single plume mask frame
%       W           = vertical velocity frame from thermOpticFlow
%
% OPTIONAL:
%       U           = horizontal velocity frame from thermOpticFlow
%       x           = x position vector
%       z           = z position vector
%       param       = parsed pulseTrack input struct - SEE
%                   'loadPulseTrackInput.m'
%
%   !!!! GO TO loadPulseTrackInput.m for all input params! !!!
%
%   For any above optional arguments left blank, they are assumed empty 
%   and given zero weight for clustering.
%
%   NAME/VALUE Pairs (entered as a struct) - see loadPulseTrackInput
%
% OUT:  
%   S = struct containing clustering and tracking results
%
%
% C Rowell March 2020

%%
    defMaxClust     = param.maxClust+1;
% narginchk(3,5)
    if nargin<7
        param = loadPulseTrackInput;
    end

    % Checks
    assert(all(size(W)==size(T)))
    assert(all(size(mask)==size(T)))
    if any(U(:))
        assert(all(size(U)==size(T)))
    end
    assert(length(x)==size(T,2))
    assert(length(z)==size(T,1))

    
    % Check pulseTracking params
    if isempty(param.lambda)
        param.lambda = defLambda;
    end
    if isempty(param.zI)
        param.zI = [1; round(sum(mask(1,:))*0.5)]; % Defaults to 1/2 mask width at base
    end
    if isempty(param.prior)
        % Default prior mask is default detection window
        [~,param.prior.mask] = getROI(mask,'iLims',[param.detectWindow(1) sum(param.detectWindow)-1],'maxRegions',1);
        
    end
    

%% Set up data array and tracking window
    % zI,detHeight,detOffset,Tprc,Gprc
    
    Tprc = param.Tpercentile;
    [X,Z] = meshgrid(x,z);
    [~,winMask] = getROI(mask,'iLims',param.zI,'jLims',param.xI,'maxRegions',1);

    % Temperature percentile cut - adjust as necessary to keep a large
    % enough sample size and ensure prior is within Temperature range
    maxPx  = sum(winMask(:)); % Max # of pixels in plume mask
    if param.minPx<maxPx
        Tcut   = prctile(T(winMask(:)),Tprc); % Temperature threshold
        keepI = and(T>=Tcut,winMask); % Pixel selection mask
        nPx    = sum(keepI(:)); % Number of pixels after cutting
        while nPx<param.minPx
            Tprc = Tprc-5; % 5 percentile increments if necessary
            Tcut   = prctile(T(winMask(:)),Tprc); % Temperature threshold
            keepI = and(T>=Tcut,winMask);
            nPx    = sum(keepI(:)); % Number of pixels after cutting
        end
    elseif param.minPx>=maxPx
        Tcut   = min(T(winMask(:)));
        Tprc   = 0;
%         keepI  = (1:maxPx)';
        keepI  = winMask;
        nPx    = maxPx;
    end
    
    % 2nd check to ensure prior is in Temperature range
    if isfield(param.prior,'Trange')
        if param.prior.Trange(1)<=Tcut
            prcrange = 100-Tprc;
            Tcut = param.prior.Trange(1);
            Tprc = (1 - sum(T(winMask)>=Tcut)./sum(winMask(:)))*100;
            Tprc_hi = min([Tprc+prcrange 100]);
            try
                Tcut_hi = max([prctile(T(winMask(:)),Tprc_hi) param.prior.Trange(2)]);
            catch
                disp('FFS')
            end  
            keepI = and(and(T>=Tcut,T<=Tcut_hi),winMask);
            nPx    = sum(keepI(:));
        end
    end
    
    % Compose cropped data vector
    imgIdx = find(keepI);
    rawDAT = [X(keepI(:)) Z(keepI(:)) T(keepI(:)) U(keepI(:)) W(keepI(:))];
    
    % ##### DAT COLUMN LABELS ######
    Labels = {'X [m]', 'Z [m]', 'T [K]', 'u [m/s]', 'w [m/s]'};
    
%% Data Selection and Standardization
  %-> Set up quick functions
    % Demean and normalize by standard deviation

    % F'n to temove nan-weighted variables, apply non-nan weights
    DatSelect = @(x,weight) x(:,~isnan(weight)).*weight(~isnan(weight));
    
    % Normalize
    S.DatMean = mean(rawDAT,1);
    S.DatStd  = std(rawDAT-S.DatMean,[],1);
    DAT = (rawDAT-S.DatMean)./S.DatStd;
    
%% Run clustering
    S.method         = 'spectral';
    S.metric         = 'euclidean';
    S.clusterWeights = param.clusterWeights;
    S.PriorWeights   = param.priorWeights;
    S.clusterMode    = param.clusterMode;
    S.Tpercentile    = Tprc;
    S.lambda         = param.lambda;
    S.nPx            = nPx;
    S.winIdx         = find(winMask);
    
    switch S.clusterMode
        case 'eig'
%             Nclust = defMaxClust;
            [~,~,S.eigd] = spectralcluster(DatSelect(DAT,S.clusterWeights),...
                defMaxClust,'Distance',S.metric);

            % Try to pick number of clusters from the eigenvalues
            Seigs = round(S.eigd,8);
            nz = Seigs~=0; % non-zero eigenvalues
            nc = 1:defMaxClust;
            belowMin = nc<param.minClust;
            aboveMax = nc>param.maxClust;
            dE = [0; diff(diff(Seigs)); 0]; % Use a crude 2nd deriv to find the biggest positive step
            dE(~nz) = 0;       % Make sure zero eigenvalues won't be chosen
            dE(belowMin) = NaN; % Zero out values below the min number of allowed clusters (max is already excluded by differencing)
            dE(aboveMax) = NaN;

            [~,S.NclustCalc] = nanmax(dE); % Biggest positive jump = estimated best number of clusters

            % Run again with estimated best number
            [S.CL,S.eigv,~] = spectralcluster(DatSelect(DAT,S.clusterWeights),...
                S.NclustCalc,'Distance',S.metric);
            
            S = getTrackedCluster(S,rawDAT,winMask,imgIdx,param); %...?
                        
        case 'fast'
            % Current plan - get chosen cluster for all allowed cluster
            % numbers (or say top 3 picks?), then optimize across those to pick the best number of
            % clusters and tracked claster (aka 2 optimizations over clusters)
            % **(remember to always exclude param.maxClust)?
            
           
            NclustRange = param.minClust:param.maxClust;

            switch length(NclustRange)>1
                case true
                    
                    Sset = S;
                    for cc = 1:length(NclustRange)
                        Sset(cc).NclustCalc = NclustRange(cc);
                        [Sset(cc).CL,Sset(cc).eigv,Sset(cc).eigd] = spectralcluster(DatSelect(DAT,Sset(cc).clusterWeights),...
                            Sset(cc).NclustCalc,'Distance',Sset(cc).metric);
                        if cc==1
                            Sset = getTrackedCluster(Sset(cc),rawDAT,winMask,imgIdx,param);
                            Sset = repmat(Sset,[length(NclustRange) 1]);
                        else
                            Sset(cc) = getTrackedCluster(Sset(cc),rawDAT,winMask,imgIdx,param);
                        end
                    end
                    Sfull    = assembleMultiCluster(Sset);

                    S = getTrackedCluster(Sfull,rawDAT,winMask,imgIdx,param);
                    S.NclustCalc = S.NclustRange(S.trackedCluster);
                    
                case false
                    
                    S.NclustCalc = NclustRange(1);
                    [S.CL,S.eigv,~] = spectralcluster(DatSelect(DAT,S.clusterWeights),...
                        S.NclustCalc,'Distance',S.metric);
                    S = getTrackedCluster(S,rawDAT,winMask,imgIdx,param); 
            end
                    
        case 'full' % formerly 'restrictive'
           %-> Basically similar to fast mode, but only allow cluster pixels
           % within uTol of old prior to count towards a new prior -
           %requires minmimizing # clusters, possibly also killing Tprc
           
             NclustRange = param.minClust:param.maxClust;

            switch length(NclustRange)>1
                case true
                    
                    Sset = S;
                    for cc = 1:length(NclustRange)
                        Sset(cc).NclustCalc = NclustRange(cc);
                        [Sset(cc).CL,Sset(cc).eigv,Sset(cc).eigd] = spectralcluster(DatSelect(DAT,Sset(cc).clusterWeights),...
                            Sset(cc).NclustCalc,'Distance',Sset(cc).metric);
                        if cc==1
                            Sset = getTrackedCluster(Sset(cc),rawDAT,winMask,imgIdx,param);
                            Sset = repmat(Sset,[length(NclustRange) 1]);
                        else
                            Sset(cc) = getTrackedCluster(Sset(cc),rawDAT,winMask,imgIdx,param);
                        end
                    end
                    Sfull    = assembleMultiCluster(Sset);

                    S = getTrackedCluster(Sfull,rawDAT,winMask,imgIdx,param);
                    S.NclustCalc = S.NclustRange(S.trackedCluster);
                    
                case false
                    
                    S.NclustCalc = NclustRange(1);
                    [S.CL,S.eigv,~] = spectralcluster(DatSelect(DAT,S.clusterWeights),...
                        S.NclustCalc,'Distance',S.metric);
                    S = getTrackedCluster(S,rawDAT,winMask,imgIdx,param); 
            end
            
            % Warp cluster
            S = restrictivePriorWarp(S,winMask,imgIdx,param);
            S.regionAvgs(S.trackedCluster,:) = mean([X(S.trackedIdx) Z(S.trackedIdx)...
                    T(S.trackedIdx) U(S.trackedIdx) W(S.trackedIdx)],1);
    end
    
    
end

function S= getTrackedCluster(S,DAT,winMask,idxList,param) %?
    

% if size(S.CL,2)==1
    % Decide on prior type (detection window or previous cluster)
    
    S.clusterAvgs = zeros(S.NclustCalc,size(DAT,2)); 
    S.regionAvgs  = zeros(S.NclustCalc,size(DAT,2)); % Largest contiguous region in each cluster
    for jj = 1:S.NclustCalc
        if size(S.CL,2)==S.NclustCalc % Full mode
            CLi = S.CL(:,jj);
        elseif size(S.CL,2)==1  % Fast mode
            CLi = S.CL==jj;
        end
        
        % Averages for each cluster
        S.clusterAvgs(jj,:) = mean(DAT(CLi,:),1);
        
        % Connected regions
        regmask = false(size(winMask));
        regmask(idxList(CLi)) = true;
        S.conncomp(jj) = bwconncomp(regmask);
        
       % Largest contiguous region for each cluster  - cluster index,
        % largest region area, region indices in input data, region averages
        [aa,kk] = max(cellfun(@(x) numel(x),S.conncomp(jj).PixelIdxList));
        S.largestClusterI(jj) = kk;
        S.clusterLargestRegion(jj) = aa;
        [~,S.regionDatIdx{jj}] = ismember(S.conncomp(jj).PixelIdxList{kk},idxList);
        S.regionAvgs(jj,:) = mean(DAT(S.regionDatIdx{jj},:),1);
        
        

    end

    % OBJECTIVE FUNCTION HERE
    objFun = getObjectiveFun(S,winMask,param);
    
%%    % ---- Sort scores for output ----
%     S.Obj.obj = objFun;
%     [S.Obj.obj, S.Obj.clusterIdx] = sort(S.Obj.obj);
%     S.Obj.Cscore = Cscore(Obj.clusterIdx);
%     S.Obj.TWDA   = Pvec(Obj.clusterIdx,:);
%     S.Obj.Pscore = Pscore(Obj.clusterIdx);
    % ---- Alt ----%
    S.trackedCluster = objFun.trackedCluster;
    [~,clusterIdx] = sort(objFun.Obj);
    S.Obj = [objFun.Cscore objFun.Pvec objFun.Pscore objFun.Obj];
    S.Obj = [clusterIdx S.Obj(clusterIdx,:)];
    % -------------    
    
    % --- TEMPORARY CODE FOR VIEWING CLUSTERS -----
%     blankmask = zeros(size(winMask));
%     T = blankmask;
%     T(idxList) = DAT(:,3);
%     blankmask = blankmask + winMask*-1;
%     blankmask = blankmask + param.prior.mask*-1;
%     for ss=1:S.NclustCalc; blankmask(S.conncomp(ss).PixelIdxList{S.largestClusterI(ss)})=ss;end % Set to cluster
% %     for ss=1:S.NclustCalc; blankmask(S.conncomp(ss).PixelIdxList{S.largestClusterI(ss)})=blankmask(S.conncomp(ss).PixelIdxList{S.largestClusterI(ss)})+1;end % Cumulative
%     [roi,~,~] = getROI(blankmask~=0);
%     [~,~,trkPol] = getROI(blankmask==S.trackedCluster,'maxRegions',1);
%     
%     
%     figure(2)
%     set(gcf,'position',[50 300 900 400]);clf
%     ax(1)=subplot(1,2,1);
%     imagesc(T)
%     axis([roi(3:4) roi(1:2)])
%     colormap(ax(1),thermgray(150))
%     set(gca,'YDir','normal')
%     hold on
%     
%     ax(2)=subplot(1,2,2);
%     imagesc(blankmask); set(gca,'YDir','normal')
%     colorbar; colormap(ax(2),jet(S.NclustCalc+3)); caxis([-2.5 S.NclustCalc+.5])
%     axis([roi(3:4) roi(1:2)])
%     hold on
%     if ~isempty(trkPol)
%         plot(ax(1),trkPol.X,trkPol.Y,'LineWidth',2)
%         plot(ax(2),trkPol.X,trkPol.Y,'LineWidth',2)   
%     end
    % ------------
    
    % Get pixel indices for tracked cluster
      %--> FILL ANY HOLES IN THE CLUSTER
    S.trackedIdx = S.conncomp(S.trackedCluster).PixelIdxList{S.largestClusterI(S.trackedCluster)};
    
end


function objFun = getObjectiveFun(S,winMask,param)

    % Distance transform for prior model
    bd = bwdist(param.prior.mask).*winMask;
    bdtop = sort(bd(:),'descend');

   % TERM 1: CURRENT CLUSTER PARAMETERS
    Cprod = prod([abs(S.regionAvgs(:,[3 5])) S.clusterLargestRegion'],2,'omitnan');
    objFun.Cscore = 1 - Cprod./max(Cprod);
    
    % TERM 2: PRIOR CLUSTER FIT
    objFun.Pvec = zeros(S.NclustCalc,4); % Add Area
    
    if isfield(param.prior,'Tsum')
        % Temperatarure measure (summed T) - incorporates Area information
        objFun.Pvec(:,1) = abs(S.regionAvgs(:,3).*S.clusterLargestRegion' - param.prior.Tsum)./param.prior.Tsum;
    elseif isfield(param.prior,'Tbar')
        % ALT (mean T):
        objFun.Pvec(:,1) = abs(S.regionAvgs(:,3)-param.prior.Tbar)'./param.prior.Tbar; %S.DatStd(3) - param.prior.Tbar;
    end
    
    if isfield(param.prior,'Wbar')    
        % Velocity measure (mean W)
        objFun.Pvec(:,2) = abs((S.regionAvgs(:,5) - param.prior.Wbar)./param.prior.Wbar);
    end
    
    if isfield(param.prior,'mask')
        % Summed distances and areas
        Dsum = zeros(S.NclustCalc,1);
        Dsum_maxes = Dsum; % Max possible summed distance for each cluster
%         clustArea = Dsum;
        for jj = 1:S.NclustCalc
            Dsum(jj) = sum(bd(S.conncomp(jj).PixelIdxList{S.largestClusterI(jj)}));
            Dsum_maxes(jj) = sum(bdtop(1:S.clusterLargestRegion(jj))) ;
        end
        
        %--- Distance normalization ----
        Dnorm = S.clusterLargestRegion'*param.pxTol; %param.uTol.*param.uMax./param.uGrid;
        % -------------------------------
        objFun.Pvec(:,3) = Dsum./Dnorm;
    end
    
    if isfield(param.prior,'Area')
        
        objFun.Pvec(:,4) = abs(S.clusterLargestRegion'-param.prior.Area)./sum(bd<param.pxTol & winMask & ~param.prior.mask,[1 2]); % Normalize to a uMax-based area change
    end

    objFun.Pscore = vecnorm(objFun.Pvec.*S.PriorWeights,2,2);
    
    % FINAL OBJECTIVE
    objFun.Obj = objFun.Cscore + S.lambda.*objFun.Pscore;
    [~,objFun.trackedCluster] = min(objFun.Obj); % Minimized objective function
        % '-> potentially allow tracking of secondary cluster when scores
        % are close? Move window with primary, fall back on secondary if
        % primary Cscore drops below certain value (ie loses heat, disappears, etc)?
        % '-> Use only in initial detection frame
end

function Sout = assembleMultiCluster(S)
    % Compile best result from multiple cluster runs for an additional
    % tracking optimization step.

    assert(and(isstruct(S),length(S)>1),'S must be struct of length>1')
    
    % Initialize universal
    Sout.DatMean        = S(1).DatMean;
    Sout.DatStd         = S(1).DatStd;
    Sout.method         = S(1).method;
    Sout.metric         = S(1).metric;
    Sout.clusterWeights = S(1).clusterWeights;
    Sout.PriorWeights   = S(1).PriorWeights;
    Sout.clusterMode    = S(1).clusterMode;
    Sout.Tpercentile    = S(1).Tpercentile;
    Sout.lambda         = S(1).lambda;
    Sout.nPx            = S(1).nPx;
    Sout.winIdx         = S(1).winIdx;
    Sout.NclustCalc     = length(S);
    Sout.NclustRange    = [S.NclustCalc];
    
    % Initialize selections
    Sout.CL                     = false(Sout.nPx,length(S));
    Sout.eigv                   = zeros(Sout.nPx,length(S));
    Sout.eigd                   = zeros(length(S),1);
    Sout.clusterAvgs            = zeros(length(S),size(S(1).clusterAvgs,2));
    Sout.regionAvgs             = zeros(length(S),size(S(1).regionAvgs,2));
    Sout.conncomp               = repmat(S(1).conncomp(1),[1 length(S)]);
    Sout.largestClusterI        = zeros(1,length(S));
    Sout.clusterLargestRegion   = zeros(1,length(S));
    Sout.regionDatIdx           = cell(1,length(S));
%     Sout.Tsum                   = zeros(1,length(S));
    
    % Loop and select best
    for ii = 1:length(S)
        Sout.CL(:,ii)               = S(ii).CL==S(ii).trackedCluster;
        Sout.eigv(:,ii)             = S(ii).eigv(:,S(ii).trackedCluster);
        Sout.eigd(ii)               = S(ii).eigd(S(ii).trackedCluster);
        Sout.clusterAvgs(ii,:)      = S(ii).clusterAvgs(S(ii).trackedCluster,:);
        Sout.regionAvgs(ii,:)       = S(ii).regionAvgs(S(ii).trackedCluster,:);
        Sout.conncomp(ii)           = S(ii).conncomp(S(ii).trackedCluster);
        Sout.largestClusterI(ii)    = S(ii).largestClusterI(S(ii).trackedCluster);
        Sout.clusterLargestRegion   = S(ii).clusterLargestRegion(S(ii).trackedCluster);
        Sout.regionDatIdx{ii}       = S(ii).regionDatIdx{S(ii).trackedCluster};
%         Sout.Tsum                   = S(ii).Tsum(S(ii).trackedCluster);
        
    end
    
end

function S = restrictivePriorWarp(S,winMask,idxList,param)

S.rawTrackedIdx = S.trackedIdx;

% Reject restrictive warp if prior has no area measure
if ~isfield(param.prior,{'Area'})
    S.offTrack = false;
    return
end

winROI = getROI(winMask);

clusterMask = false(size(winMask));
clusterMask(S.trackedIdx) = true;

priorMask = param.prior.mask(winROI(1):winROI(2),winROI(3):winROI(4));
clustMask = clusterMask(winROI(1):winROI(2),winROI(3):winROI(4));

outwardCut = (bwdist(priorMask)>param.pxTol).*clustMask;
inwardCut  = and(bwdist(~priorMask)>param.pxTol,~clustMask);
clustMask  = clustMask + inwardCut - outwardCut;
% Fill any holes
clustMask  = imfill(clustMask);

% Quick smoother to eliminate small offshoots
sN = 3;
clustMask  = convn(double(clustMask),ones(sN,sN)./sN.^2,'same')>0.5;

% Where we wander off track, just employ a stopping condition for now. Can
% consider implemnting a re-iteration approach in future.
if ~any(clustMask(:))
    warning('Tracked cluster has no overlap with prior.')
    S.offTrack = true;
    return
else
    S.offTrack = false;
end

clusterMask(winROI(1):winROI(2),winROI(3):winROI(4)) = clustMask;

% Update fields
S.trackedIdx = find(clusterMask);
[~,S.regionDatIdx{S.trackedCluster}] = ismember(find(clusterMask),idxList);
S.clusterLargestRegion(S.trackedCluster) = length(S.trackedIdx);

end
        
