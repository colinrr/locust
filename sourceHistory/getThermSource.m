function S = getThermSource(D,V,varargin) %trackParams,detParams,plotflag)
%  Function to scan a static region of thermal image over time, retrieving
%  statistics and detecting pulses of hot material.
% IN:
%   D = thermal data cube struct from "getThermCube.m"
%   V = velocity data cube struct from "thermOpticFlow.m"
%
% NAME/VALUE arguments (can be entered as struct with corresponding fieldnames):
%   trackParams = struct controlling the window region with FIELDS:
%       trackWindowStart  :  bottom pixel of main window from which thermal 
%                        statistics are retrieved. Default = 1;     
%       trackWindowHeight : vertical extent (pix) of main window 
%                           Default = equivalent of 50 m
%       removeAtmo        : true/false flag. (default false). Where true
%                           and atmo struct is present in D, atmospheric 
%                           profile will be removed. Plume masks will be
%                           modified to remove pixels below a threshold of
%                           -nHalfMax*(atmo.Tmax-atmo.T_halfMax),
%                           where nHalfMax is typically 2. "nHalfMax" can be
%                           set manually as an additional Name/Value pair.
%
%      det_params = secondary struct of STA/LTA detection parameters array, DEFAULT = [5 40 2 1.6 0]
%         FIELDS:  (DEFAULTS: 0.05, 20, 2, 20, 2, 1.6, 0, 3, 'continuous',  respectively)
%         --> taperLength - length (fraction of input signal) of taper 
%                       applied to both end of signals (Tukey Window)
%         --> movMinLength - Length (SECONDS) of moving minima window to
%               for signal pre-processing signal. This can
%               have a major influence on the number/period of detected
%               signals and should be chosen with care. A good value is 
%               typically slightly larger than the dominant oscillation 
%               period of the signal.
%         --> l_sta    - STA window length (s)
%         --> l_lta    - LTA window length (s)
%         --> th_on    - STA/LTA trigger on threshold
%         --> th_off   - STA/LTA trigger off threshold
%         --> min_dur  - Minimum event duration to be recorded (s)
%         --> lta_mode - (string) ... default: 'continuous'
%            '-> 'frozen' - LTA window fixEVENT_ON && ((sta_to_lta(count) <= th_off) || count == length(y)) || no_signal(count)ed in place after trigger is turned on 
%                    while STA window continues forward.
%            '-> 'continuous' - LTA window continues w/ STA window after trigger 
%                        is turned on (Same behavior as before trigger)
%
%       outputFile      : Path to save source detection data
%
% A SEPARATE DETECTION WINDOW IS DEPRECATED - the 2 parameters below
% can be safely ignored and will be removed.
%       detWindowHeight   : height (pixels) of detection window
%       detWndowOffset    : offset height (pixels) of detection window
%   --> Because velocities and detections are at the leading edge of a moving
% pulse, but thermal statistics should be retrieved from the pulse body,
% the detection window should be vertically offset above the main window.
% Default detection window is the upper half of the main window:
% (ie detWindowOffset  = round(trackWinHeight * 0.5) - tracking starts half a window height up) 
% (ie detWindowScale  = round(trackWinHeight * 0.5)  - tracking window height is one half the stats window) 
%
% OUT:
%   S = struct containing statistics for static 'source' window (typcally
%       plume region just above the vent). With fields:
%       T0: Thermal stats output (from getThermStats) for source window
%       V0: Velocity stats output (from getThermStats) for source window
%       Td:
%       Vd:
%       tTrig: [nx2] start, stop times of n detected sources
%       iTrig: [nx2] start, stop indices of n detected sources
%       trackParams: same as input, populated with defaults where necessary
%
% C Rowell March 2020


fprintf('\n========= Thermal Detection and Source Statistics =========\n')

%%
if nargin<2
    V = [];
end

defTrackWinStart    = 1;
defTrackWinHeight   = round(100/D.dx);
% defDetWinOffset     = 0; %round(defTrackWinHeight/2); % (pixels)
% defDetWinHeight     = defTrackWinHeight; %round(defTrackWinHeight/2); % (pixels)
defDetParams = struct('l_sta',[], 'l_lta',[], 'th_on', [], 'th_off', [],...
    'min_dur',[],'lta_mode',[],'plotflag',[]); 
defDetChan = 5; % Default channel to use for detection in 'prctile' data type
defDataType = 'var';
defPrcVals  = [5 25 50 75 95];

p = inputParser;
addParameter(p,'trackWindowStart',defTrackWinStart)
addParameter(p,'trackWindowHeight',defTrackWinHeight)
addParameter(p,'removeAtmo',false)
% addParameter(p,'detectionWindowOffset',defDetWinOffset)
% addParameter(p,'detectionWindowHeight',defDetWinHeight)
addParameter(p,'detParams',defDetParams)
addParameter(p,'dataType',defDataType)
addParameter(p,'detectionChannel',defDetChan)
addParameter(p,'prcvals',defPrcVals)
addParameter(p,'satVal',[])
% addParameter(p,'plotflag',0)

% The following 2 are a holdover from deprecated code, but are still used as
% defaults and could be re-implemented
addParameter(p,'windowDuration',1)
addParameter(p,'nFramesOverlap',0)
addParameter(p,'nHalfMax',2)
addParameter(p,'outputFile',[])
parse(p,varargin{:});

S = p.Results;

    %% Source function - thermal max, variance, mean/rms/median?

    M = size(D.T,1);
    N = size(D.T,2);
    P = size(D.T,3);

    dt = mean(diff(D.t)); % Time spacing
    
    % REMOVE ATMO PROFILE?
%     if and(S.removeAtmo,isfield(D,'atmo'))
%         [D.T,D.mask] = removeAtmoProfile(D.mask,D.T,D.atmo,S.nHalfMax,true);
%         D.T = D.T - D.atmo.Tinterp - D.atmo.T_halfMax;
%         newMask = D.T >= -S.nHalfMax*(D.atmo.Tmode-D.atmo.T_halfMax);
%         D.mask = and(D.mask,newMask);
%     end

    % Get time windows
    [src_tI,~] = getSTFTColumns(P,S.windowDuration,S.nFramesOverlap,1/dt); % Something weird with this wintTime

    Pwin             = size(src_tI,2); % Number of time windows

    S.trackWindow = [S.trackWindowStart S.trackWindowHeight];
%     S.detectWindow = [S.trackWindowStart+S.detectionWindowOffset S.detectionWindowHeight];
    src_zI = repmat((S.trackWindow(1):sum(S.trackWindow)-1)', [1 Pwin]); % Set source region
%     det_zI = repmat((S.detectWindow(1):sum(S.detectWindow)-1)', [1 Pwin]); % Set detection window
%     S = rmfield(S,{'trackWindowStart','trackWindowHeight','detectionWindowOffset','detectionWindowHeight'});
    S = rmfield(S,{'trackWindowStart','trackWindowHeight'});
    
    % Get source window stats
    if and(S.removeAtmo,isfield(D,'atmo'))
        S.T0 = getThermStats(D.T,D.mask,D.z,D.t,src_tI,src_zI, S.prcvals, D.atmo, S.satVal);
    else
        S.T0 = getThermStats(D.T,D.mask,D.z,D.t,src_tI,src_zI, S.prcvals, [], S.satVal);
    end
    if ~isempty(V)
        S.V0 = getThermStats(V.Vz,D.mask,V.z,V.t,src_tI,src_zI, S.prcvals);
    else
        S.V0 = [];
    end
    % [V0f,V0spec] = filterVelocities(V0,dt,'low',3,4,true);

    % Detection window stats
%     S.Td = getThermStats(D.T,D.mask,D.z,D.t,src_tI,det_zI, S.prcvals);
%     S.Vd = getThermStats(V.Vz,D.mask,D.z,V.t,src_tI,det_zI, S.prcvals);
%     [Vdf,Vdspec] = filterVelocities(Vd,dt,'low',2,4,filter_plotflag);
    % Tmax = squeeze(max(D.T.*D.mask,[],2));


    %% Implement detection here
    disp('Running source pulse detection...')
    
    switch S.dataType
        case 'prctile'
            Tdetect = S.T0.(S.dataType)(S.detectionChannel,:);
        otherwise
            Tdetect = S.T0.(S.dataType);
    end
    [S.tTrig,S.iTrig,S.detParams] = thermPulseDetection(Tdetect, S.T0.t, S.detParams);

    %% Data save
    if ~isempty(S.outputFile)
        if exist(S.outputFile,'file')
            warning(sprintf('Output file already exists:\n\t%s',S.outputFile))
            cmd = input('Overwrite? (y/n) ','s');
            if strcmp(cmd,'y')
                disp('Overwriting...')
                writeout = true;
            else
                disp('Abort save')
                writeout = false;
            end
        else
            writeout = true;
        end
        
        if writeout
            fprintf('Writing source stats:\n\t%s\n',S.outputFile)
            save(S.outputFile,'S')
        end
    end
    
end