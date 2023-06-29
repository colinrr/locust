function param = loadPulseTrackInput(varargin)
% Parse pulseTrack Name/Value paris into struct of parameter settings.
% Parameters can be entered as a struct with corresponding field names, or
% as standard Name/Value argument pairs.
%
% OPTIONS:
%
%   'trackWindow'     : [I H] INITIAL vertical position of TRACKING window.
%                     --> I =  bottom pixel position
%                     --> H =  height (in pixels) of tracking window.
%                     Alternate: [Nx2] matrix, where N = length(iTrigger), to use a
%                     unique tracking window for each track
%
%   'detectWindow'    : same as param.trackWindow, but for the DETECTION window if
%                     thermPulseDetection was used. This provides an initial prior 
%                     model for cluster detection.
%                     --> Default = upper half of tracking window.
%
%   'iTrigger'        : starting frame indices for each track (ie detection frames)
%                     --> scalar or vector of indices along 3rd DIM of Vz
%                     --> this is an output of thermPulseDetection.m
%
%   'Tpercentile'     : Percentile value with which to filter out cool 
%                     temperatures during clustering and tracking. 
%                     --> Default = 70.
%
%   'lambda'          : Regularization parameter for cluster optimization
%
%   'maxClust'        : Maximum number of clusters allowed
%
%   'minClust'        : Minimum number of clusters allowed
%
%   'clusterWeights'  : Variable weights for clustering....DEPRECATED
%                       -->           [X   Z   T    U   W]
%                       --> Default = [0.5 0.5 1.5  1   2];
%   'PriorWeights'    : Variable weights for OBJECTIVE PRIOR TERM.
%                       -->           [Temp  W_vel  Dist  Area]
%                       --> Default = [0.5   0.25   2     0.5];
%
%   'clusterMode'     : 'full', 'fast', or 'eig'. 
%                       --> 'full' uses in iterative technique to optimize 
%                       output over the range of possible cluster sizes 
%                       (minClust:maxClust).
%                       --> 'fast' automatically calculates the optimum
%                       number of clusters based on the relative sizing of
%                       tracking window and tracked features
%                       --> 'eig' estimates the based number of clusters
%                       based on the spectral clustering eigenvalues. It is
%                       generally only used on the first time step and is
%                       not intended to be specified manually. It is
%                       included as an option for testing purposes.
%
%   'prior'           : "Prior model" - struct of parameters from of 
%                       previous cluster, typically left empty for input,
%                       and initial mask is generated from 'detectWindow'.
%                       Any values input here will be used as the INITIAL
%                       prior for the first frame.
%                          Nclust : Number of clusters used in prior
%                          mask   : prior cluster after velocity warp.
%                                   Used to calculate distance transform
%                          Tbar   : Average temperature of prior cluster
%                          Wbar   : Average vertical velocity of prior cluster
%                          Area   : "area" (# pixels) in prior cluster
%   'maskTracking'    : true [default] or false. If true, when the initial
%                       detection window contains the plume mask top, the
%                       tracking window will automatically follow the plume
%                       mask top, keeping the window at a size similar to
%                       the plume mask (assumes plume front is feature of
%                       interest). If false, tracks as normal in all cases.
%
%   OTHER PARAMS (can be tuned manually with input, but not typically
%                 necessary as they are all calculated internally):
%   'nHalfMax'   : After removing atmospheric profile, a data threshold will be set 
%                  at -nHalfMax*(Tmax-T_halfMax), where (Tmax-T_halfMax) is
%                  the half-width of the temperature distribution for
%                  temperature values that fit the atmospheric profile.
%                  (see fitAtmoProfile to make sense of this).
%                  Pixels below this value will be discarded as background. 
%   'minPx'      : Minimum number of pixels to include in clustering.
%                  Tpercentile will be dynamically adjusted when the amount
%                  of input data is too small.
%   'winSzRatio' : Target ratio of vertical cluster extent to window
%                  extent. Initially set by tracking window height by
%                  default.
%   'uMax'       :  A typical maximum velocity for any given frame. Normally
%                   calculated internally.
%   'uGrid'      : Characteristic grid speed, dx/dt, calculated internally.
%   'uTol'       :  Tolerance multiplier for maximum tracking velocity (used
%                   in Prior calculation). uMax and uTol together set a
%                   hard limit for the velocity of tracking window and soft
%                   limit for the movement of tracked clusters.
%   'pxTol'      : Maximum motion distance (units of pixels) corresponding
%                  to uTol. Calculated internally, but can be forced here.
%   'memoryN'    : Number of previous frames over which to build prior.
%
% C Rowell Jun 2020
 

%% Default tracking params
    trackWindow    = []; % zI?
    detectWindow   = [];
    iTrigger       = 1;
    Tpercentile    = 30;
%     zI             = [];
    clusterMode    = 'full';
    maxClust       = 6;
    minClust       = 3;
    clusterWeights = [0.5 0.5 2   1   1.5];
    priorWeights   = [0.5 0.25 2 0.5];
    prior          = [];
    maskTracking   = true;
    lambda         = 0.4;
    nHalfMax       = 2;
    minPx          = 256;
    winSzRatio     = 1.5;
    uMax           = [];
    uGrid          = [];
    uTol           = 2.5;
    pxTol          = [];
    memoryN        = [];
    stopTime       = [];
    stopHeight     = [];
    applyAtmo      = true;
    qcplot         = false;
    
%% Parse input
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true;
    
    addParameter(p,'trackWindow',trackWindow)
    addParameter(p,'detectWindow',detectWindow)
    addParameter(p,'iTrigger',iTrigger)
    addParameter(p,'Tpercentile',Tpercentile)
%     addParameter(p,'zI',zI)
    addParameter(p,'clusterMode',clusterMode)
    addParameter(p,'maxClust',maxClust)
    addParameter(p,'minClust',minClust)
    addParameter(p,'clusterWeights',clusterWeights)
    addParameter(p,'priorWeights',priorWeights)
    addParameter(p,'prior',prior)
    addParameter(p,'maskTracking',maskTracking)    
    addParameter(p,'lambda',lambda)
    addParameter(p,'nHalfMax',nHalfMax)
    addParameter(p,'minPx',minPx)
    addParameter(p,'winSzRatio',winSzRatio)
    addParameter(p,'uMax',uMax)
    addParameter(p,'uGrid',uGrid)    
    addParameter(p,'uTol',uTol)
    addParameter(p,'pxTol',pxTol)
    addParameter(p,'memoryN',memoryN)
    addParameter(p,'stopTime',stopTime)
    addParameter(p,'stopHeight',stopHeight)  
    addParameter(p,'applyAtmo',applyAtmo)
    addParameter(p,'qcPlot',qcplot)
%     addParameter(p,'',)
    parse(p,varargin{:})
    
    param = p.Results;
    puf = fieldnames(p.Unmatched);
    if ~isempty(puf)
        warning('Unrecognized pulseTrack Name/Value argument(s):')
        fprintf('\t%s\n',puf{:})
    end

end