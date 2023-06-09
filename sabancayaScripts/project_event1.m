
% =======================================================================
%                        pulseTracker Project
% =======================================================================
% C Rowell, March 2020

% ANALYSIS FOR EVENT 1 (18-05-25A, sequence 4)
%(1-180525AA-04: May 25, Continuous plume, 08:01-08:06)

disp('Initializing Event 25A-4...')
% clear all; close all
clearvars -except D V; close all

addpath(genpath('plumeTracking'))
addpath(genpath('pulseTracker'))
addpath(genpath('preProcTools'))
addpath(genpath('ijcv_flow_code'))

% =========================== USER INPUT =============================
%                             vvvvvvvvvv
%% ================== DATA DIRECTORIES ==================

homeDir = '.'; % LINUX/MAC
dataDir   = fullfile(homeDir,'testData/');

% -------- MAIN WORKING directory for this event --------
thermDir   = fullfile(dataDir,'event1/');

% -------- Directory of RAW .MAT files for individual frames --------
matDir  = fullfile(thermDir,'mat/');
headsIRB = fullfile(matDir,'params-heads.mat');
headsIRT = fullfile(matDir,'params.mat');
matHeads = fullfile(matDir,'frameHeadsFixed.mat');

% -------- FIGURE OUTPUT --------
figdir = fullfile(dataDir,'Saba_manuscript_data/mat-figs/');

allIdx = 1:2439; % All frame indices - useful reference/plot tool
%% ======== Image Registration (thermPixelReg) ==========
regDir  = fullfile(thermDir,'reg-mat/');
regHeads = fullfile(regDir,'frameHeads.mat');
% regParams = fullfile(regDir,'registration_params_2020-03-17_n1929.mat');
regParams = fullfile(regDir,'registration_params_2020-03-23_n1929.mat');
% regParams = fullfile(regDir,'registration_params_2020-03-15_n1929.mat');

regMode = 'full';

regROI = [1 1024 692 768]; % [x1 x2 y1 y2] - Use this to register a select image region (ie foreground)
refIdx = 510; % Index of reference frame (other frames will be transformed to this one)

% Subset of frames that require registration
% regIdx = [511:769]; 
regIdx = 511:2439;

 % Temperature range for scaling registration images. Final images are not
 % scaled, but registration needs simple grayscale image
% regTscale = [200 419.85]; % Kelvin
regTscale = [260 350]; % Kelvin

%% ====== INITIAL PLUME CALCULATIONS (mapPixels and plumeTrackCalcs) ======
% !!!!! CAREFUL !!!!!
% USE ORIGINAL IMAGES AND DIMENSIONS ONLY FOR INITIAL mapPixels run
% -> NOT! REGISTERED IMAGES. 
% IF REGISTERING, OBTAIN refPix FROM THE SAME IMAGE USED FOR refIdx.

satVal = 420; % Saturation brightness temp for this imagery
            %'-> Not certain, no saturated values observed, IRB calibration dependent.

% Distance mapping
%           lat         lon         elev (m)
obsLLE = [-15.738265, -71.84283, 5243]; % Camera coords using DEM elevation

% vent   = [-15.786744, -71.855919, 5911]; % Best guess from SRTM? or GE?
vent   = [-15.787498, -71.856122, 5919]; % New guess from alos12m DEM. Based on lowest crater point along LOS to BI explosion first jet
% vent_utm = [1.93989e+05 8.252493e+06 5.919e+03];
hfov = 32.36; vfov=24.55;


% Control point to calibrate camera azimuth/tilt
%           lat         lon         elev (m)
% ref_utm = [1.943832634403001e+05 8.252621640823159e+06 5954];
ref_utm = [194376 8252622 5954];
% refLLE = [-15.740833, -71.852500, 5913.000000]; % Landmark coords
[refLLE(1),refLLE(2)] = utm2deg(ref_utm(1),ref_utm(2),'19 L'); refLLE(3) = ref_utm(3);

% Use registered images for these references IF registration was done
% refPix = [681 295]; % [y x] pixel coords in ref image corresponding to refLLE
refPix = [682 296];
% imsz = [764 1020]; % Registered image size (or raw if no registration was done)
imsz = [768 1024]; % Raw images
% !!!!!!!!!!!!!!!!!!!!

% DEM plots/calcs
demfile = fullfile(dataDir,'dem_alos/alos12m_saba_clip60x60_utmZ19.tif');
% dem_roi = {[192500 198000],[8251500 8258000]};
dem_roi = {[193500 194500],[8252000 8253000]};
% dem_roi = {[185000 205000],[8245000 8264000]};

plotCalcs = true;
plot_image_projection = true;

geomf  = fullfile(thermDir,'geometry_2020-07-19_registered.mat');
%% ============= THERMAL PRE-PROCESSING (preProcThermal) =============
% created pre-processes thermal images to make plumeTracker's job easier
 % - design foreground masks, and/or time-varying smooth heaviside threshold
 % function to scale out background/noise

procDir = fullfile(thermDir,'preProc-mat/'); % pre-processed frames for plumeTracker
procIdx = 1:30;
% ---- FOREGROUND MASKING ----
thresh_fg   = 262; % Temperature threshold applied to REFERENCE image to cut out foreground
nullVal     = 210;    % Min value to scale temperatures down to

% Polygon design
polys = false; % false = EXCLUDE pixels from FOREGROUND region
verticalDisplace = 6; % Manually shift foreground boundary up by this many pixels as a buffer

fgFile  = fullfile(procDir,'foreground_mask');
fgPolyFile = fullfile(procDir,'foreground_poly.mat');

% --- HEAVISIDE FILTERING -----
ITW = [1 ; ... % Frame Index for control points.
       230 ; ... % T0 - Center of smooth threshold filter, Kelvin
       10 ];        
% ITW = []; 
HtFile  = fullfile(procDir,'Temp_time_hist_fgmask.mat');

% SWITCHES
get_foreground_mask = false; % Design a single mask file to remove foreground
view_histograms     = true;  % Useful for designing filters
% test_filters   = 
apply_filters         = false; % Apply "filters" to frames
artificial_ref        = false; % Artificially zero out reference image to nullVal
 
%% ============= PLUME TRACKER (mainTrackPlume) ==============
% INPUT: most recent data set of: matDir, regDir, or procDir
% -------- Plume Track OUTPUT directory --------
ptInputDir = procDir;
ptInputHeads = regHeads;
plumeTrackDir = fullfile(thermDir,'PTresults/');

            
ref = 1; % Index of reference image to use
deb = 1; % Starting image? Need not be the same as reference image
          % >> USE VECTOR to ignore fin and dN, and use indices in vector
          % instead
          
fin = 30; % Leave empty to grab all files to end
% deb = 500;fin=670; 
% freq=10;   % Sample frequency, Hz
dN = 1;    % Track every dN frames
deb = [deb:dN:allIdx(end)];
% deb = [deb:dN:fin];
procIdx = [ref deb];
fixpix = [413 687]; % Fixed pixel at the plume base center [x y]

% Plotting
delT = [220 350]; % Temperature range in Kelvin to plot vids
% delT = [0 3]; % For gradient images
plume_params = true; % Plots height/width tracker bars in the images
png   = false;
gif   = false;
video = true;

ptHeads = fullfile(plumeTrackDir,'plumeTrack_output.mat'); % plumeTrack Output params
polyFile = [];

%% =========== THERMAL IMAGE REGRIDDING (interpThermal) ===========
interpDir = fullfile(thermDir,'interp-mat-Jul2020/');
% interpDir = fullfile(thermDir,'interp-mat-test/');
interpIn = regDir;

interpIdx = [];

vmax = 100; % Maximum velocity allowed in plume mask movement
interpHeads = fullfile(interpDir,'frameHeads.mat'); % -> interpThermal ouput

%% ============== THERMAL DATA CUBE SETUP ================
cubeDir = fullfile(thermDir,'thermCubeAnalysis/');
% thermIdx = [];
thermIdx = 1:2:2439; % Reduced frames to expirement with opticalflow
% thermIdx = 1:5:2439; % Reduced frames (EVEN MORE) to expirement with opticalflow

%  [x1 x2 y1 y2] % Must be based off INTERPOLATED IMAGES
ROI = [315 938 1 704]; % REALLY big window to catch the whole plume, let masking do the work of cropping

% Current Therm data cube file
thermCube   = fullfile(cubeDir,'thermStats_2020-07-20_z704_x624_t2439.mat');
% thermCube   = fullfile(cubeDir,'thermStats_2020-10-27_z704_x624_t1220.mat'); % Fewer frames to try with optical flow

%% === ATMOSPHERIC PROFILE CHARACTERIZATION/REMOVAL (get/fitAtmoProfile) ===
% MODIS/AIRS files
atmoHDF = {
            fullfile(dataDir,'MODIS_atmospheric_profiles/MYD07_L2.A2018145.1750.061.2018146174008.hdf');
            fullfile(dataDir,'MODIS_atmospheric_profiles/MOD07_L2.A2018145.1505.061.2018146014305.hdf');
            fullfile(dataDir,'AIRS_atmospheric_profiles/AIRS.2018.05.24.188.L2.RetStd_IR.v6.0.31.1.G19090065901.hdf');
%             fullfile(dataDir,'AIRS_atmospheric_profiles/AIRS.2018.05.25.056.L2.RetStd_IR.v6.0.31.1.G19090102252.hdf');
%             fullfile(dataDir,'AIRS_atmospheric_profiles/AIRS.2018.05.25.056.L2.RetSup_IR.v6.0.31.1.G19090102252.hdf')
%             fullfile(dataDir,'AIRS_atmospheric_profiles/AIRS.2018.05.26.186.L2.RetStd_IR.v6.0.31.1.G19090134728.hdf')
            };
        
atmoRefFrames = 10:25:2350;
% atmoRefzI     = [];
atmoImgROI     = [196 478 1 333]; % Cut out values below 500m - too much variation from vent effects
localTimezoneOffset = -5;
atmoTthresh    = 240;

% Large portion of plume masks are outside of main plume, so plume scale
% and off-axis position (and therefore Z error) will be largely
% overestimated. So filter data primarily using the tight ROI above and
% allow zError Uncertainty threshold to be large.
TerrorThresh       = 1.5;
zUncertaintyThresh = 700;

% atmoFile = fullfile(dataDir,'MODIS_atmospheric_profiles/atmo_25B.mat');
atmoFile = fullfile(thermDir,'atmoProfiles_x1_20-10-27.mat');
atmoRefCube = thermCube; % Which thermal data to use for atmospheric reference
%% ============== DATA CUBE VELOCIMETRY (thermOpticFlow2O) ================

opticData = thermCube;

opticParams.Sub = []; % 3rd-dimension subscript (time) in data cube to select frames for optic flow.
opticParams.maskPad = 64; % Def=15; Algorithm will crop to field around mask, plus a pixel pad this thick
opticParams.method = [];  % 'method' input for 'estimate_flow_interface'
opticParams.flowParams = {}; % Additional input name-value pairs for 'estimate_flow_interface'

% Filtering
% lowpass cutoff period (seconds)
OFfiltlim = 1.25;  % For ALL frames (step of 1)
OFfiltlim = 2.5;  % For every 5th frame (step of 5)OFfiltlim = 1.25; % lowpass cutoff period (seconds)
filtTest = true;
filtPlot = true;

% Name out output velocity file
% velCube = fullfile(cubeDir, 'opticFlowCNL_20-07-31_n2439.mat'); % Unfiltered
velCube = fullfile(cubeDir, 'opticFlowCNL_20-07-31_n2439_filtered.mat'); % Filtered
%% ============== THERMAL SOURCE DETECTION (getThermSource) ===============

muWinHeight = 30; % Height of windows for mean image profile

% STA/LTA data input
thermSource.dataType      = 'var'; % Which field of T0 (source window thermal stats) to use for detection
% thermSource.dataType      = 'prctile'; 
thermSource.detectionChannel = 5;
% thermSource.prcvals          = [5 20 50 80 95];
thermSource.prcvals          = [5 25 50 75 95];
thermSource.satVal           = satVal;
thermSource.removeAtmo       = true;

% Source and detection window parameters
thermSource.trackWindowStart = 1;
thermSource.trackWindowHeight = 20;
% thermSource.detectionWindowOffset = 0; % (pixels)
% thermSource.detectionWindowHeight = 20; % (pixels)

% STA/LTA detection parameters
detParams.taperLength   = [];
% detParams.movMinLength  = 15;
detParams.movMinLength  = 10;
detParams.l_sta         = [];
detParams.l_lta         = [];
detParams.th_on         = 1.0;
detParams.th_off        = 0.7; 
detParams.min_dur       = []; 
detParams.lta_mode      = [];
detParams.plotflag      = true;
thermSource.detParams = detParams;

% detection_plotflag  = true; % Plot detection results
% filter_plotflag     = true; % Plot pre- and post-filter velocity w/ spectra

% Source detection file
% -> source stats and final detection triggers (from
% manual_tracking_25A4)
% thermSourceFile = fullfile(cubeDir,'thermSource_25A4_2020-09-30.mat');
thermSourceFile = fullfile(cubeDir,'thermSource_25A4_2020-10-29.mat');
%% ============== THERMAL PULSE TRACKING (pulseTrack) ================

% rawTrackFile = fullfile(cubeDir,'25A4_tracks_raw.mat'); % Track output
rawTrackFile = fullfile(cubeDir,'Tracks_25A4_combined.mat'); % Track output
trackFile = fullfile(cubeDir,'Tracks_25A4_processed.mat'); % Track output

allTrackFile = fullfile(dataDir,'allTracks_24A_25A4_25B.mat');

trackExcelFile = fullfile(dataDir,'Data_and_Event_index.xlsx');
xlsSheet       = 'pulseTrack_Notes';

%% ============== Time Averaged Images (getAveragedImage) ================
avgROI = [1 272 1 268];

% Frame indices to use for generating averaged image
avgIdx = []; % Full history

avgWinHeight = 30;
% avgWinOver   = 29;

% avgImg = fullfile(cubeDir,'time_averaged_all_T_Vx_Vz.mat');
% avgImg = fullfile(cubeDir,'tAvg_25A4_steady_2021-01-10.mat');
% avgImg = fullfile(cubeDir,'tAvg_25A4_steady_2022-07-07.mat');
avgImg = fullfile(cubeDir,'tAvg_25A4_steady_2022-07-16.mat');

%% RUN pulseTracker
% 1) Pre-process thermal
% preProc25A

% 2) 
% pulseTrackDriver

% 3) Run feature tracking
% manual_tracking_25A4
