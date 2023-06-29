    
% =======================================================================
%                   pulseTracker Project, Event 2
% =======================================================================
% C Rowell, July 2018

% ANALYSIS FOR EVENT 2  (18-05-24A) 
%(1-180524AA-00: May 24, 10:30 explosion)
% clear all; close all;
clearvars -except D V; %close all
% =========================== USER INPUT =============================
disp('Initializing Event 24A...')

addpath(genpath('plumeTracking'))
addpath(genpath('pulseTracker'))
addpath(genpath('preProcTools'))
addpath(genpath('ijcv_flow_code'))
%% ================== DATA DIRECTORIES ==================

setHomeDir

% -------- MAIN WORKING directory for this event --------
thermDir   = fullfile(dataDir,'event2/');

% -------- Directory of RAW .MAT files for individual frames --------
% matDir  = fullfile(thermdir,'mat/');
matDir  = fullfile(thermDir,'mat/');
matHeads = fullfile(matDir,'params.mat');


procDir = fullfile(thermDir,'mat/grad-scale2/'); % pre-processed frames for plumeTracker
outputDir = fullfile(procDir,'PTresults/');

asc_params = fullfile(thermDir,'RT_1030.txt');
% asc_params = [];

% glob_spec = '*.tif';
glob_spec = 'RT_1030_corrected_*.txt';


% -------- FIGURE OUTPUT --------
figdir = fullfile(thermDir,'figures');

%% ======== Image Registration (plumePixelReg) ========

regParams = []; % No Registration done
%% ============= THERMAL PRE-PROCESSING (preProcThermal) =============
% created pre-processes thermal images to make plumeTracker's job easier
 % - design foreground masks, and time-varying smooth heaviside T threshold
 % function to scale out background/noise
  %  - > Separate parameter key for this right now

%% ====== INITIAL PLUME CALCULATIONS (MapPixels and plumeTrackCalcs) ======
% !!!!! CAREFUL !!!!!
% USE ORIGINAL IMAGES AND DIMENSIONS ONLY FOR INITIAL mapPixles run
% -> NOT! REGISTERED IMAGES. 
% IF REGISTERING, OBTAIN refPix FROM THE SAME IMAGE USED FOR refIdx.


% Distance mapping
%           lat         lon         elev (m)
obsLLE = [-15.738265, -71.84283, 5243]; % Camera coords using DEM elevation
vent   = [-15.787498, -71.856122, 5919]; % New guess from alos12m DEM. Based on lowest crater point along LOS to BI explosion first jet
hfov = 32.36; vfov=24.55;

% Control point to calibrate camera azimuth/tilt
ref_utm = [194376 8252622 5954];
[refLLE(1),refLLE(2)] = utm2deg(ref_utm(1),ref_utm(2),'19 L'); refLLE(3) = ref_utm(3);
refPix = [712 347]; % [y x] pixel coords in ref image corresponding to refLLE
imsz = [768 1024]; % Raw images

% DEM plots/calcs
demfile = fullfile(dataDir,'dem_alos/alos12m_saba_clip60x60_utmZ19.tif');
% dem_roi = {[192500 198000],[8251500 8258000]};
dem_roi = {[193500 194500],[8252000 8253000]};

plotCalcs = true;
plot_image_projection = false;

% Camera/target positioning geometry
geomf  = fullfile(thermDir,'geometry_2020-06-26.mat');

%% ============= PLUME TRACKER (mainTrackPlume) ==============
satVal = 424.7; % Saturation brightness temp for this imagery

ref = 370; % Index of reference image to use
deb = 382; % Starting image? Need not be the same as reference image
          % >> USE VECTOR to ignore fin and dN, and use indices in vector
          % instead -> STILL IMPLEMENTING
fin = 520;
% fin = 1576; % Leave empty to grab all files to end
% freq=10;   % Sample frequency, Hz
dN = 2;    % Track every dN frames
deb = [378:2:1573];
fixpix = [465 717]; % Fixed pixel for plume base [x y]

% Plotting
delT = [230 400]; % Temperature range in K to plot vids
% delT = [0 3]; % For gradient images
plume_params = true; % Plots height/width tracker bars in the images
png   = false;
gif   = false;
video = true;

ptHeads = fullfile(matDir,'plumeTrack_output.mat'); % plumeTrack Output params
polyFile = fullfile(thermDir,'manual_polygons_all.mat');

%% =========== THERMAL IMAGE REGRIDDING (interpThermal) ===========
interpIn = matDir;
interpIdx = [378:1572];

vmax = 100;

% interpDir = fullfile(thermDir,'interp-mat/');
% interpDir = fullfile(thermDir,'interp-mat2/');
interpDir = fullfile(thermDir,'interp-mat-Jun2020/');
interpHeads = fullfile(interpDir,'frameHeads.mat'); % -> interpThermal ouput
%% ============== THERMAL DATA CUBE SETUP (getThermCube) ================
% Data cubes for thermal and velocity data
cubeDir = fullfile(thermDir,'thermCubeAnalysis/');

thermIdx = [378:1572]';
% thermIdx = [378:395]';

%  [x1 x2 y1 y2] % Must be based off INTERPOLATED IMAGES
 % Big window to catch the whole plume, let masking do the work of cropping
ROI = [340 917 80 723];

% Current Therm data cube file
thermCube   = fullfile(cubeDir,'thermStats_2022-07-07_z644_x578_t1195_revCorr.mat'); % Fixed atmo profile
%% === ATMOSPHERIC PROFILE CHARACTERIZATION/REMOVAL (get/fitAtmoProfile) ===
% MODIS/AIRS files
atmoHDF = {
            fullfile(dataDir,'AIRS_atmospheric_profiles/AIRS.2018.05.24.066.L2.RetStd_IR.v6.0.31.1.G19090062535.hdf')
            };
        
atmoRefFrames = (750:5:1110);
atmoImgROI = [90 644 1 578]; % Crop out an area of images for best fit
atmoTthresh = 240;           % Use a temperature threshold to assist mask clipping
localTimezoneOffset = -5;

TerrorThresh       = 1.0;
zUncertaintyThresh = 1000; % (T error thresh is the more important limiter)

atmoFile = fullfile(thermDir,'atmoProfiles_x2_20-10-23.mat');
atmoRefCube = thermCube; % Which thermal data to use for atmospheric reference
%% ============== DATA CUBE VELOCIMETRY (thermOpticFlow2O) ================
% -----------OPTICAL FLOW------------
opticData = thermCube;

opticParams.Sub = []; % 3rd-dimension subscript (time) in data cube to select frames for optic flow.
opticParams.method = [];
opticParams.maskPad = 64; % Def=15; Algorithm will crop to field around

% Filtering
OFfiltlim = 2; % lowpass cutoff period (seconds)
filtTest = false;
filtPlot = true;

% Name out output velocity file
velCube = fullfile(cubeDir, 'opticFlowCNL_20-06-29_n1195_filtered.mat'); % Filtered
%% ============== THERMAL SOURCE DETECTION (getThermSource) ===============

muWinHeight = 30; % Height of windows for mean image profile

% STA/LTA data input
% thermSource.dataType      = 'var'; % Which field of T0 (source window thermal stats) to use for detection
thermSource.dataType      = 'prctile'; 
thermSource.detectionChannel = 5;
% thermSource.prcvals          = [5 20 50 80 95];
thermSource.prcvals         = [5 25 50 75 95];
thermSource.satVal          = satVal;
thermSource.removeAtmo      = true;

% Source and detection window parameters
thermSource.trackWindowStart = 1;
thermSource.trackWindowHeight = 20;

% STA/LTA detection parameters
detParams.taperLength   = [];
% detParams.movMinLength  = 15;
detParams.movMinLength  = 10;
detParams.l_sta         = [];
detParams.l_lta         = [];
detParams.th_on         = 0.7;
detParams.th_off        = []; 
detParams.min_dur       = []; 
detParams.lta_mode      = [];
detParams.plotflag      = true;
thermSource.detParams = detParams;

detection_plotflag  = true; % Plot detection results
filter_plotflag     = true; % Plot pre- and post-filter velocity w/ spectra

% Source detection file
% -> source stats and final detection triggers (from
% manual_tracking_25A4)
thermSourceFile = fullfile(cubeDir,'thermSource_24A_2020-10-29.mat');

%% ============== THERMAL PULSE TRACKING (pulseTrack) ================

% Structure tracking inputs are contained in the separate file:
% sabancayaScripts/structureTracking_event2

% Output track files
rawTrackFile = fullfile(cubeDir,'Tracks_24A_raw2.mat');
trackFile    = fullfile(cubeDir,'Tracks_24A_processed.mat');

allTrackFile = fullfile(dataDir,'allTracks_24A_25A4_25B.mat');

%% ============== Time Averaged Images (getAveragedImage) ================
% avgROI = [1 400 46 366];
avgROI = [1 290 46 250];

% Frame indices to use for generating averaged image
% avgIdx = 1:2:686; % Full history, no decay regime (~ 0 - 150 s)
avgIdx = 216:686; % Quasi-steady regime (~47-150 s)

avgWinHeight = 15;
avgWinOver   = 14;

avgImg = fullfile(cubeDir,'tAvg_24A_quasiSteady_2022-07-16.mat');
%% ------- Rise Diagrams -------
ventPix     = [473 716];  % [x y] plume base (there must always be a mask pixel at this y value)
profileType = 'max';
Ridx        = [380:2:1573]'; %[15:2:41]'; 

%% RUN pulseTracker
% pulseTrackDriver

% manual_tracking_24A