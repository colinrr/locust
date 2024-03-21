% =======================================================================
%                        pulseTracker Project
% =======================================================================
% C Rowell, March 2020

% ANALYSIS FOR EVENT 25A, sequence 0
%(1-180525AA-04: May 25, Continuous plume, 07:35:18-07:035:14 local)

disp('Initializing Event 25A-0...')
% clear all; close all
% clearvars -except D V; close all

addpath(genpath('plumeTracking'))
addpath(genpath('pulseTracker'))
addpath(genpath('preProcTools'))
addpath(genpath('ijcv_flow_code'))

% =========================== USER INPUT =============================
%                             vvvvvvvvvv
%% ================== DATA DIRECTORIES ==================

homeDir = 'C:\Users\crowell\';  % SlimJim
% homeDir = 'D:\';    % FB
% homeDir = '~'; % LINUX/MAC
dataDir   = fullfile(homeDir,'/Kahuna/data/sabancaya_5_2018/');

% -------- MAIN WORKING directory for this event --------
thermDir   = fullfile(dataDir,'image_exports/25A-0/');

% -------- Directory of RAW .MAT files for individual frames --------
matDir  = fullfile(thermDir,'mat/');
IRBheads = fullfile(matDir,'frameHeads_crop_to_dropped_frames.mat'); % Good timestamps - best guess was that 3 frames were dropped so N=1139
IRTheads = fullfile(matDir,'frameHeads_IRT_bad_timestamps.mat'); % Bad timestamps attached to data, N=1136
matHeads = fullfile(matDir,'frameHeadsFixed.mat'); % OFFICIAL frame headers here

% -------- FIGURE OUTPUT --------
figdir = fullfile(thermDir,'figures/');

allIdx = 1:2439; % All frame indices - useful reference/plot tool

%% ======== Image Registration (thermPixelReg) ==========
regROI = [1 1024 692 768]; % [x1 x2 y1 y2] - Use this to register a select image region (ie foreground)

regParams = [];
% NO REGISTRATION PERFORMED

%% ====== PROJECTION CALCULATIONS (mapPixels and plumeTrackCalcs) ======

satVal = 420; % Saturation brightness temp for this imagery
            %'-> Not certain, no saturated values observed, IRB calibration dependent.

% Distance mapping
%           lat         lon         elev (m)
obsLLE = [-15.738265, -71.84283, 5243]; % Camera coords using DEM elevation
vent   = [-15.787498, -71.856122, 5919]; % guess from alos12m DEM. Based on lowest crater point along LOS to BI explosion first jet

% Camera field of view (degrees)
hfov = 32.36; vfov=24.55;

% Control point to calibrate camera azimuth/tilt
%           lat         lon         elev (m)
ref_utm = [194376 8252622 5954];
[refLLE(1),refLLE(2)] = utm2deg(ref_utm(1),ref_utm(2),'19 L'); refLLE(3) = ref_utm(3);

% Use registered images for these references IF registration was done
refPix = [685 299];
imsz = [768 1024]; % Raw images
% !!!!!!!!!!!!!!!!!!!!

plotCalcs = true;
plot_image_projection = true;

geomf  = fullfile(thermDir,'geometry_2024-03-18.mat');

%% ============= THERMAL PRE-PROCESSING (preProcThermal) =============
% created pre-processes thermal images to make plumeTracker's job easier
 % - design foreground masks, and/or time-varying smooth heaviside threshold
 % function to scale out background/noise
 % ==> SEE preProcWorkflow_25A0.m for workflow steps

procDir = fullfile(thermDir,'preProc-mat/'); % pre-processed frames for plumeTracker
heads = matHeads; % fullfile(matDir,'frameHeads.mat');

fgFile  = fullfile(procDir,'foreground_mask.mat');

ref_idx     = 1; % Reference frame to use
procIdxTest = 1:30; % Initial set to test plumeTracker on with filters
procIdx     = [];   % Full set

thresh_fg   = 261; % Temperature threshold applied to REFERENCE image to cut out foreground
satVal      = 420; % Saturation temperature for thermal vid
nullVal     = 210;    % Min value to scale temperatures down to

% % --- HEAVISIDE FILTERING -----
ITW = [1 ; ... % Frame Index for control points.
       230 ; ... % T0 - Center of smooth threshold filter, Kelvin
       10 ]; 

% ---- PLOT SWITCH
plot_flag = false;
 
 %% ============= PLUME TRACKER (mainTrackPlume) ==============
% INPUT: most recent data set of: matDir, regDir, or procDir
% -------- Plume Track OUTPUT directory --------
ptInputDir = procDir;
ptInputHeads = matHeads;
plumeTrackDir = fullfile(thermDir,'PTresults/');

            
ref = 1; % Index of reference image to use
deb = 1; % Starting image? Need not be the same as reference image
          % >> USE VECTOR to ignore fin and dN, and use indices in vector
          % instead
          
fin = 1136; % Leave empty to grab all files to end
% deb = 500;fin=670; 
% freq=10;   % Sample frequency, Hz
dN = 1;    % Track every dN frames
deb = [deb:dN:allIdx(end)];
% deb = [deb:dN:fin];
procIdx = [ref deb];
fixpix = [424 690]; % Fixed pixel at the plume base center [x y]

% Plotting
delT = [220 350]; % Temperature range in Kelvin to plot vids
% delT = [0 3]; % For gradient images
plume_params = true; % Plots height/width tracker bars in the images
png   = false;
gif   = false;
video = true;
PTplotFlags = [plume_params png gif video];

ptHeads = fullfile(plumeTrackDir,'plumeTrack_output.mat'); % plumeTrack Output params
polyFile = [];

%% =========== THERMAL IMAGE REGRIDDING (interpThermal) ===========
interpDir = fullfile(thermDir,'interp-mat/');
interpIn = matDir;

interpIdx = [];
interpIdxTest = 1:30;

vmax = 60; % Maximum velocity allowed in plume mask movement
interpHeads = fullfile(interpDir,'frameHeads.mat'); % -> interpThermal output

%% ============== THERMAL DATA CUBE SETUP ================
cubeDir = fullfile(thermDir,'thermCubeAnalysis/');

thermIdx = []; % All frames
thermIdxAtmo = round(linspace(1,1136,201)); % To make a subset for atmo profile fitting

%  [x1 x2 y1 y2] % Must be based off INTERPOLATED IMAGES
ROI = [325 942 1 709]; % REALLY big window to catch the whole plume, let masking do the work of cropping

% Current Therm data cube file
thermCube   = fullfile(cubeDir,'thermStats_2024-03-19_z709_x618_t1136.mat');

%% === ATMOSPHERIC PROFILE CHARACTERIZATION/REMOVAL (get/fitAtmoProfile) ===
% MODIS/AIRS files
atmoHDF = {
            fullfile(dataDir,'MODIS_atmospheric_profiles/MYD07_L2.A2018145.1750.061.2018146174008.hdf');
            };
        
atmoRefFrames = []; %10:25:2350;

% atmoRefzI     = [];
atmoImgROI     = [124 530 1 365]; % Cut out values below 500m - too much variation from vent effects
localTimezoneOffset = -5;
atmoTthresh    = 240;

% Large portion of plume masks are outside of main plume, so plume scale
% and off-axis position (and therefore Z error) will be largely
% overestimated. So filter data primarily using the tight ROI above and
% allow zError Uncertainty threshold to be large.
TerrorThresh       = 3; %1.5;
zUncertaintyThresh = 1000; %700;

atmoFile = fullfile(thermDir,'atmoProfiles_x1_24-03-20.mat');
atmoRefCube = fullfile(cubeDir,'thermStats_2024-03-20_z709_x618_t201.mat'); % Which thermal data to use for atmospheric reference

%% ============== THERMAL SOURCE DETECTION (getThermSource) ===============

muWinHeight = 30; % Height of windows for mean image profile

% STA/LTA data input
thermSource.dataType      = 'var'; % Which field of T0 (source window thermal stats) to use for detection
% thermSource.dataType      = 'prctile'; 
thermSource.detectionChannel = 5;
thermSource.prcvals          = [5 25 50 75 95];
thermSource.satVal           = satVal;
thermSource.removeAtmo       = true;

% Source and detection window parameters
thermSource.trackWindowStart = 1;
thermSource.trackWindowHeight = 20;

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
thermSourceFile = fullfile(cubeDir,'thermSource_25A0_2024-03-20.mat');