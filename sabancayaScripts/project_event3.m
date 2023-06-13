
% =======================================================================
%                        pulseTracker Project
% =======================================================================
% C Rowell, March 2020

% ANALYSIS FOR EVENT 1 (18-05-25B)
% (BI052500: May 25, 15:11 'big bang')
% =========================== USER INPUT =============================
% clear all; close all
clearvars -except D V;  %close all
disp('Project: Event 25B...')

%% ================== DATA DIRECTORIES ==================

setHomeDir

% -------- MAIN WORKING directory for this event --------
thermDir   = fullfile(dataDir,'/event3/');

% -------- Directory of RAW .MAT files for individual frames --------
matDir  = fullfile(thermDir,'mat-Mar20/');
headsIRB = fullfile(matDir,'params-heads.mat');
headsIRT = fullfile(matDir,'params.mat');
matHeads = fullfile(matDir,'frameHeadsFixed.mat');

% -------- FIGURE OUTPUT --------
figdir = fullfile(dataDir,'figures/');

%% ======== (2) Image Registration (thermPixelReg) ==========
regDir  = fullfile(thermDir,'reg-mat/'); % Directory to registered frames
regHeads = fullfile(regDir,'frameHeads.mat'); % Frame header information
regParams = fullfile(regDir,'registration_params_2020-03-23_n2575.mat');

regMode = 'full';

regROI = [1 1024 705 768]; % [x1 x2 y1 y2] - Use this to register a select image region (ie foreground)
refIdx = 7; % Index of reference frame (other frames will be transformed to this one)

% Subset of frames that require registration
regIdx = [1:6 8:2576];

 % Temperature range for scaling registration images. Final images are not
 % scaled, but registration needs simple grayscale image
regTscale = [230 330]; % Kelvin

%% ====== INITIAL PLUME CALCULATIONS (mapPixels and plumeTrackCalcs) ======
% !!!!! CAREFUL !!!!!
% USE ORIGINAL IMAGES AND DIMENSIONS ONLY FOR INITIAL mapPixels run
% -> NOT! REGISTERED IMAGES. 
% IF REGISTERING, OBTAIN refPix FROM THE SAME IMAGE USED FOR refIdx.

satVal = 403.1; %403.15; % Saturation brightness temp for this imagery

% Distance mapping
%           lat         lon         elev (m)
obsLLE = [-15.750128, -71.836807, 5168]; % Camera coords, but using DEM elevation

vent   = [-15.787498, -71.856122, 5919]; % New guess from alos12m DEM. Based on lowest crater point along LOS to BI explosion first jet
hfov = 32.36; vfov=24.55;

% Reference point to calibrate camera azimuth/tilt
ref_utm = [194376 8252622 5954]; 
[refLLE(1),refLLE(2)] = utm2deg(ref_utm(1),ref_utm(2),'19 L'); refLLE(3) = ref_utm(3);
refPix = [692 362]; % new guess from /mat-Mar20/25B_tau1_em1_0007.mat
imsz = [768 1024]; % Raw image dimension
% !!!!!!!!!!!!!!!!!!!!

plotCalcs = true;
plot_image_projection = false;

% Camera/target positioning geometry
geomf  = fullfile(thermDir,'geometry_2020-06-02_registered.mat'); % Corrected projection, adjusted for registration

%% ============= PLUME TRACKER (mainTrackPlume) ==============
% INPUT: most recent data set of: matDir, regDir, or procDir
% -------- Plume Track OUTPUT directory --------
plumeTrackDir = fullfile(thermDir,'PTresults/');

ref = 1; % Index of reference image to use
deb = 7; % Starting image? Need not be the same as reference image
          % >> USE VECTOR to ignore fin and dN, and use indices in vector
          % instead
fin = 2578; % Leave empty to grab all files to end
dN  = 2;    % Track every dN frames

deb = [(deb:10)';(11:dN:1121)' ;(1123:4:fin)'];
fixpix = [517 698]; % New for registered images

% Plotting
delT = [230 400]; % Temperature range in K to plot vids
plume_params = false; % Plots height/width tracker bars in the images
png   = false;
gif   = false;
video = true;

PTplotflags = [plume_params,png,gif,video];

ptHeads = fullfile(regDir,'plumeTrack_output.mat'); % plumeTracker Output params

polyFile = [];

%% =========== THERMAL IMAGE REGRIDDING (interpThermal) ===========
interpIdx = 11:1300; % all
interpIdx = 11:110; % demo subset
interpIn = regDir; % Image source directory

vmax = 100; % Maximum velocity allowed in plume mask movement

interpDir = fullfile(thermDir,'interp-mat/');
interpHeads = fullfile(interpDir,'frameHeads.mat');

%% ============== THERMAL DATA CUBE SETUP (getThermCube) ================
cubeDir = fullfile(thermDir,'thermCube/');
% thermIdx = (11:1299)';
thermIdx = (11:110)'; % demo cube

% Image region-of-interest to crop
ROI = [270 840 1 710]; % ALL plume up to Idx 1300

% Last output therm data cube file
thermCube   = fullfile(cubeDir,'thermStats_2020-10-26_z710_x571_t1289.mat'); % This and previous were functionally identical, just overwrote 1st by accident


% thermFile   = fullfile(cubeDir,'thermStats_2020-03-25_z710_x571_t1294.mat'); % OLD FILE

% thermCube   = fullfile(cubeDir,'thermStats_2020-06-08_z710_x571_t1289.mat');
thermCube   = fullfile(cubeDir,'thermStats_2020-10-26_z710_x571_t1289.mat'); % This and previous were functionally identical, just overwrote 1st by accident

% Atmo profile reference cube
refCube  = fullfile(cubeDir,'thermStats_2020-06-08_z710_x571_t343.mat');
% refCube  = fullfile(cubeDir,'thermStats_2020-05-22_z710_x571_t29.mat');
% Long time series ref cube
atmoRefCube  = fullfile(cubeDir,'thermStats_2020-06-09_z710_x571_t205.mat');
slimCube     = fullfile(cubeDir,'thermStats_2020-06-20_z710_x571_t205.mat');
% atmoRefCube  = fullfile(cubeDir,'thermStats_2020-06-20_z710_x571_t430.mat');

%% === ATMOSPHERIC PROFILE CHARACTERIZATION/REMOVAL (get/fitAtmoProfile) ===
% MODIS/AIRS files
atmoHDF = {
            fullfile(dataDir,'MODIS_atmospheric_profiles/MYD07_L2.A2018145.1750.061.2018146174008.hdf');
%             fullfile(dataDir,'MODIS_atmospheric_profiles/MOD07_L2.A2018145.1505.061.2018146014305.hdf');
            fullfile(dataDir,'AIRS_atmospheric_profiles/AIRS.2018.05.25.056.L2.RetStd_IR.v6.0.31.1.G19090102252.hdf');
%             fullfile(dataDir,'AIRS_atmospheric_profiles/AIRS.2018.05.25.056.L2.RetSup_IR.v6.0.31.1.G19090102252.hdf')
%             fullfile(dataDir,'AIRS_atmospheric_profiles/AIRS.2018.05.26.186.L2.RetStd_IR.v6.0.31.1.G19090134728.hdf')
            };
        
% atmoRefFrames = 994:50:1294;
atmoRefFrames = 167:205;
% atmoRefzI     = [];
% atmoRefzI     = 140:710; % Cut out values below 500m - too much variation from vent effects
atmoTthresh     = 235;           % Use a temperature threshold to assist mask clipping
atmoImgROI      = [140 710 1 571];
localTimezoneOffset = -5;

TerrorThresh       = 1.0;
zUncertaintyThresh = 350;

% atmoFile = fullfile(dataDir,'MODIS_atmospheric_profiles/atmo_25B.mat');
atmoFile = fullfile(thermDir,'atmoProfiles_x2_20-10-23.mat'); % Updated w/ fit statistics

%% ============== DATA CUBE VELOCIMETRY (thermOpticFlow2O) ================
% -----------OPTICAL FLOW------------
% Choose input data
opticPlotFrames = false;
opticData = thermCube;

% opticFlow parameters
opticParams.Sub = []; % 3rd-dimension subscript (time) in data cube to select frames for optic flow.
% opticParams.Sub = 290:320; % Small set for OpticFlow2O test
opticParams.maskPad = 64; % Def=15; Algorithm will crop to field around mask, plus a pixel pad this thick
opticParams.method = [];  % 'method' input for 'estimate_flow_interface'
opticParams.flowParams = {}; % Additional input name-value pairs for 'estimate_flow_interface'

% Filtering
OFfiltlim = []; %1.2; % lowpass cutoff period (seconds)

% Name out output velocity file
% velCube = fullfile(cubeDir, 'opticFlowCNL_20-04-09_n1294.mat'); % OLD FILE
% velCube = fullfile(cubeDir, 'opticFlowCNL_20-06-18_n1289.mat'); % Raw
velCube = fullfile(cubeDir, 'opticFlowCNL_20-06-18_n1289_filtered.mat'); % Filtered
%% ============== THERMAL SOURCE DETECTION (getThermSource) ===============
% STA/LTA data input
thermSource.dataType        = 'var'; % Which field of T0 (source window thermal stats) to use for detection
% thermSource.dataType      = 'prctile'; 
thermSource.detectionChannel = 4;
% thermSource.prcvals          = [5 20 50 80 95];
thermSource.prcvals         = [5 25 50 75 95];
thermSource.satVal          = satVal;
thermSource.removeAtmo      = true;

% Source and detection window parameters
% thermSource.trackWindowStart = 25;
% thermSource.trackWindowHeight = 45;
% thermSource.detectionWindowOffset = 15; % (pixels)
% thermSource.detectionWindowHeight = 25; % (pixels)
thermSource.trackWindowStart = 20;
thermSource.trackWindowHeight = 25;
% thermSource.detectionWindowOffset = 0; % (pixels)
% thermSource.detectionWindowHeight = 25; % (pixels)


% STA/LTA detection parameters
detParams.taperLength   = [];
% detParams.movMinLength  = 15;
detParams.movMinLength  = 10;
detParams.l_sta         = [];
detParams.l_lta         = [];
detParams.th_on         = 0.44;
detParams.th_off        = 0.44; 
detParams.min_dur       = []; 
detParams.lta_mode      = [];
detParams.plotflag      = true;
thermSource.detParams = detParams;

% Source detection file
% -> source stats and final detection triggers (from
% manual_tracking_25A4)
% thermSourceFile = fullfile(cubeDir,'thermSource_25B_2020-10-01.mat');
thermSourceFile = fullfile(cubeDir,'thermSource_25B_2020-10-29.mat');
%% ============== THERMAL PULSE TRACKING (pulseTrack) ================

% NAME/VALUE Params for manual tracking
% TrigI = 2; % Which detections from getThermSource to track (optional)
% trackPar.Tpercentile = 40;
% trackPar.lambda = 0.4;

% trackPar.iTrig = iTrig(2:4,1);
% trackPar.iTrigger       = 160;
% trackPar.iTrigger       = 140;
% trackPar.trackWindow    = [thermSource.trackWindowStart thermSource.trackWindowHeight];
% trackPar.detectWindow   = [trackPar.trackWindow(1)+thermSource.detectionWindowOffset thermSource.detectionWindowHeight];
% trackPar.detectWindow   = [16 25];
% trackPar.trackWindow    = trackPar.detectWindow;
% trackPar.minClust       = 4;
% trackPar.maxClust       = 6;
trackPar.uMax           = 26.82215; % More efficient
trackPar.memoryN        = 9;        % More efficient
% trackPar.clusterMode    = 'restrictive';

% Specified for testing purposes
% trackPar.clusterWeights = [1 1 2   1   1.5];
% trackPar.priorWeights   = [0.5 0.25 2 0.5];
% trackPar.uTol           = 2.5;

% Pulse Tracking parameters
% trackParams.detectionWindowOffset = 25; % (pixels)
% trackParams.detectionWindowHeight = 15; % (pixels)
% trackParams.trackWindowStart = 1;
% trackParams.trackWindowHeight = 40;
% trackParams.Tpercentile             = 70;
% trackParams.Gpercentile              = 90;

% detection_plotflag  = true; % Plot detection results
% filter_plotflag     = true; % Plot pre- and post-filter velocity w/ spectra

% Output track files
rawTrackFile = fullfile(cubeDir,'v2Tracks/Vtracks_25B_raw.mat');
trackFile    = fullfile(cubeDir,'v2Tracks/Vtracks_25B_processed.mat');

allTrackFile = fullfile(dataDir,'allTracks_24A_25A4_25B.mat');

%% ============== Time Averaged Images (getAveragedImage) ================
% avgROI = [1 450 1 392];
% avgROI = [1 312 1 392]; % First 1000 meters
% avgROI = [1 243 1 392]; % First 800 meters
avgROI = [];

% Frame indices to use for generating averaged image
% avgIdx = 244:840; % First clearance of 1000 m onward to almost fully decayed
% avgIdx = 155:840; % First clearance of 800 m onward to almost fully decayed

% Main emissions stage (~0-26.7 s)
avgIdx = 1:265; 
avgROI = [1 450 1 392];
avgImg = fullfile(cubeDir,'tAvg_25B_early_t0-27_2022-11-02.mat');

% --- Decaying regime () (26.7-84.7 s)
% avgIdx = 320:2:840; 
% avgROI = [1 450 1 392];
% % avgImg = fullfile(cubeDir,'tAvg_25B_decay_t27-85_2022-07-07.mat');
% avgImg = fullfile(cubeDir,'tAvg_25B_decay_t27-85_2022-07-16.mat');

% --- Full history (~ 0 - 84.7 s)
% avgIdx = 1:2:840; 
% avgROI = [1 450 1 392];
% % avgImg = fullfile(cubeDir,'tAvg_25B_all_t0-87_2022-07-07.mat');
% avgImg = fullfile(cubeDir,'tAvg_25B_all_t0-87_2022-07-16.mat');

% --- All significant pulses (~ 0 - 55 s)
% avgIdx = 1:2:540; 
% avgROI = [1 450 1 392];
% avgImg = fullfile(cubeDir,'tAvg_25B_allLargePulses_2022-07-07.mat');
% avgImg = fullfile(cubeDir,'tAvg_25B_allLargePulses_2022-07-16.mat');


avgWinHeight = 20; % Used for smoothing retrieved height tracks
% avgWinOver   = 20; % No longer used

% avgImg = fullfile(cubeDir,'tAvg_25B_initialPulse.mat');
% avgImg = fullfile(cubeDir,'tAvg_25B_allLargePulses.mat');
% avgImg = fullfile(cubeDir,'tAvg_25B_allLargePulses_2021-01-10.mat');

%% RUN pulseTracker
% pulseTrackDriver

% pickVelocityTesting
% run sandbox_scripts/pulseTrack_ensemble
% manual_detections_25B

% vidGenerator