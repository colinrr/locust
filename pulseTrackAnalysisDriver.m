% =======================================================================
%                pulseTracker post-Processing and analysis
% =======================================================================
% C Rowell, Jan 2021

% Workflow and driver for processing of pulseTracker combined dataset
clear all; close all
%% Build path
addpath(genpath('.'))
%%

setHomeDir
trackDir  = fullfile(dataDir,'pulseTrack_analysis/');
figDir    = fullfile(dataDir,'figures/');

% Velocity and thermal cubes for `buildTrackDataset.m'. Requires figShare
% data repository
thermCubes = {
    fullfile(dataDir,'figShare/thermalDataCubes/thermStats_2020-06-29_z644_x578_t1195_revCorr.mat');
    fullfile(dataDir,'figShare/thermalDataCubes/thermStats_2020-07-20_z704_x624_t2439.mat');
    fullfile(dataDir,'figShare/thermalDataCubes/thermStats_2020-10-26_z710_x571_t1289.mat');
             };
         
velCubes = {
    fullfile(dataDir,'figShare/opticalFlowDataCubes/opticFlowCNL_20-06-29_n1195_filtered.mat');
    fullfile(dataDir,'figShare/opticalFlowDataCubes/opticFlowCNL_20-07-31_n2439_filtered.mat');
    fullfile(dataDir,'figShare/opticalFlowDataCubes/opticFlowCNL_20-06-18_n1289_filtered.mat');
             };
         
% For creating mask polgyons, not strictly needed
maskFile = fullfile(dataDir,'pulseTrack_analysis/allMasks_3events.mat');  % Plume masks only

% RAW TRACK FILE output from 'buildTrackDataset.m'
tkFile           = fullfile(trackDir,'allTracks_raw_24A_25A4_25B_2022-07-20.mat');

%% File and track selection for scaleTracks (R, T curve fitting)

% Time averaged images
avgImgFile   = fullfile(dataDir,'pulseTrack_analysis/eventAveragedImages_2022-07-16.mat');

% Which tracks?
trackSet = []; % Get em all
trackSet2 = [1 20]; % Leading fronts of Event 2 and 3 respectively in raw track set

% FOR SELECTION OF RADIUS AND TEMPERATURE POINTS FOR CURVE FITTING,
% see pulseTrack_processingV3.xlsx.  
%  -->These numbers are derived from manually examining each track for 
% quality in terms of how well it follows a feature, how visible/complete 
% that feature is, how distorted that feature becomes in wind, etc. 
% Indices therefore encompass regions where tracking quality is GOOD and 
% captures the evolution of an ENTIRE structure


% Z ERROR THRESHOLDS FOR IMAGE FILTERING (zErrThresh, zErrUncThresh).
% Default dTzero for now.
zErrThresh = [300 500]; %[150 300]; % A relaxed constraint (allows more data, not much qualtiy cost)
zErrThresh2 = [300 500]; % Specific to quality data in leading fronts of Event 2 and 3

% getSensitivity = true;

%%  UPDATED scaleTracks INPUT FIELDS
trackScalingFile = fullfile(dataDir,'pulseTrack_analysis/tracks_gauss_polys_n26_2022-08-04-B.mat');
% testSet = [1 2 13 20 22 24];
testSet = 20; %[1 9 20 25];

Rcurve              = 'npxInt';
Tstruct             = 'trackStatsInt';
Tcurve              = 'prctile';
filtFractionThresh  = 0.9;
satFractionThresh   = 0.1;
getSensitivity      = false;

% Track 20 test values
RzLim = [539 1175]; 
TzLim = RzLim;
Ridx  = [];
Tidx  = [];

%% Track file updated with...

% Reconstructed images, gauss fits, and profile retrievals - intermediate
% steps that are not included in demo or figShare
trackDataFile = fullfile(trackDir,'trackData_26tracks_2022-12-05.mat');
trackPolyFile = fullfile(dataDir,'pulseTrack_analysis/tracks_gauss_polys_n26_2022-12-05.mat');

% Averaged image curve fitting results
avgImgFile = fullfile(dataDir,'pulseTrack_analysis/AverageImageScalingResults_V3.mat');


%% scaleTracks input runs
trackUnsteadyFile1 = fullfile(dataDir,'pulseTrack_analysis/trackFitResults_unsteady_v1_22-12-17.mat');

% Common input for curve fitting
sIc.filtFracThresh  = 0.2;
sIc.satFracThresh   = 0.1;
sIc.dT0Limits       = [-10 10];% Absolute c error limits about 0 (K) . Tighter test: [-5 2]; 
sIc.cLimits         = [-0.2 0.2]; % Relative c error limits about 0

cLimitsRange = [26:-2:0]'.*[-1 1]; % To test leading front sensitivity, Events 2/3

% Default T measure = 95th prctile by default

% Average z0 values - an improvement?? Nope.
load(avgImgFile,'AvgImgT')
z0Avg = [repmat(AvgImgT.z0_combo(2,:),[8 1]); repmat(AvgImgT.z0_combo(1,:),[11 1]); repmat(AvgImgT.z0_combo(3,:),[7 1])];
% z0Avg = [];
z0_tkAvg = [repmat([-150.2221  -55.4562   39.3096],[8 1])
            repmat([-225.8152 -109.2311    7.3530],[11,1])
            repmat([-41.6822   -2.1720   37.3382],[7 1])];
        
% Iteration set; , using ensemble average track z0        
load(trackUnsteadyFile1,'z0adjusted')
z0adjusted = z0adjusted([12:19 1:11 20:26],:); % Have to reverse the revised order here
% original track number indices in new event order

% ---- Curve Fitting Runs -----

% Re-arranged to capture all - THIS IS THE FINAL RESULT FOR MANUSCRIPT
sI1.desc      = 'All-main';
sI1.Rcurve    = 'npxInt';     % Choice of radius measure
sI1.Tstruct   = 'trackImg';   % Choice of temperature dataset
sI1.z0        = [];
sI1.Tcurve    = 'prctile';
sI1.trackSet  = [9:19 1:8 20:26]; %[2:19 21 23:26]; % test set
sI1.oFile     = 'trackFitResults_v2_22-12-06'; % MAIN RESULTS SET, SABA MS

% For a few cases where appropriate, used z0 estimate combining information
% from mask fit - THIS IS THE IMPROVED FIT TO THESE FOUR TRACKS THAT ARE USED
sI2.desc      = 'main-subset';
sI2.Rcurve    = 'trkCombined';     % Choice of radius measure
sI2.Tstruct   = 'trackImg';   % Choice of temperature dataset
sI2.z0        = [];
sI2.Tcurve    = 'prctile';
sI2.trackSet  = [1 10 20 21]; %[2:19 21 23:26]; % test set
sI2.oFile     = 'trackFitResults_R-trkCombined';
tkNoOut = [12 2 20 21]; 

% Sensitivity analysis on z0
sI3.desc      = 'z0sensitivity';
sI3.Rcurve    = 'npxInt';     % Choice of radius measure
sI3.Tstruct   = 'trackImg';   % Choice of temperature dataset
sI3.z0        = [];
sI3.Tcurve    = 'prctile';
sI3.trackSet  = [1 3 12 15]; %[2:19 21 23:26]; % test set
sI3.sensitivityFlag = true;
sI3.oFile     = 'trackFitResults_z0sens_22-12-17'; % Use with sI3 - testing B sensitivity to z0

% Re-arranged to capture all
sI4.desc      = 'z0_unsteady_adjusted';
sI4.Rcurve    = 'npxInt';     % Choice of radius measure
sI4.Tstruct   = 'trackImg';   % Choice of temperature dataset
sI4.z0        = z0adjusted;
sI4.Tcurve    = 'prctile';
sI4.trackSet  = [9:19 1:8 20:26]; %[2:19 21 23:26]; % test set
sI4.oFile     = 'trackFitResults_fixedz0_22-12-17'; % Use with sI4 - fixed z0 using ensemble track average 

% z0 from time-averaged images 
sI5.desc      = 'z0_tAvg';
sI5.Rcurve    = 'npxInt';
sI5.Tstruct   = 'trackImg';
sI5.z0        = z0Avg; %(sI10.trackSet,:);
sI5.Tcurve    = 'prctile';
sI5.trackSet  = [9:19 1:8 20:26];
sI5.oFile     = 'trackFitResults_fixedz0-tAvg_23-06-13'; % Use with sI5 - fixed time-averaged image z0

% Which of the above runs to do?
sIuse = [sI1 sI2]; % main MS results
% sIuse = [sI3]; % z0 sensitivity tests
% sIuse = [sI4]; % Testing track ensemble average z0 - generally good results in agreement with main set, just wider scatter
% sIuse = [sI5]; % Testing time-averaged z0 - extremely negative and probably unphysical results
% sIuse = [sI1 sI4 sI5]; % comparing

avgErrMode = 2; % Method to calculate average error across tracks
% 1 = mean and std.dev. of central estimates
% 2 = weighted mean using rmse results

% Output figure names
scaleFig = 'f09_trackResults_v3';
zBfig    = 'sf_ZvsB_v3';
suppTkFig = 'sf_trackFits_v2';

%% WORKFLOW Switches

flag_gatherAveragedImages   = false; % Gather time-averaged images into a combined dataset
flag_scaleAveragedImages    = false; % Get scaling results for gathered time-averaged images
flag_combineTrackDataset    = false; % Get raw tracks from all events together in 1 dataset
flag_buildTrackDataset      = false; % Track post-processing: reconstructed images, gaussian fits, statistics profiles 
flag_getTrackPolygons       = false; % Calculate track outline polygons - total and frame-by-frame
flag_scaleTracks            = true; % Perform curve fitting calculations to get drdz, z0, T^B decay exponent
    save_scalingOutput          = false;
flag_unsteady               = false;

% Figures
plotResults  = true;
print_scalingResults = false;
print_sampleTrackResults = false;
%% DO THE THING - WORKFLOW

% (1) ------ R,Texp time-averaged images for track comparisons ------ 
if flag_gatherAveragedImages
    gatherAveragedImages
end
if flag_scaleAveragedImages
    scaleAverageImages % -> Needs dr/dz, condensed output, and saving.
end

% (2) ------ Combine tracking datasets from tracker output for different events ------
if flag_combineTrackDataset
    combineTrackDataset 
end

% (3) ------ Get reconstructed "images" and Z profiles of t-dependent tracks ------
if flag_buildTrackDataset
    [tk,allTrkPar] = buildTrackDataset(tkFile,thermCubes,velCubes,trackSet,zErrThresh(1,:),winLength,trackDir); % Get full set
    % Specific processing for Event 2 and 3 starting pulses
%     [tk,allTrkPar] = buildTrackDataset(tkFile,thermCubes,velCubes,trackSet2,zErrThresh2(1,:),winLength,trackDir); % Get full set
end

% (4) ------ Get track outline polygons ------ (these are mainly for visualization and QC - eg. Supplementary tracking figure)
if flag_getTrackPolygons
    tk = getTrackPolygons(trackDataFile,maskFile,trackPolyFile,false);
end

% (5) ------ ANALYSIS: Run retrieval of scales, spreading rate, and exponent ------
if flag_scaleTracks
    load(trackPolyFile,'tk')
    
    nS = length(sIuse);
    
    % Looping over scaling input sets
    for si = 1:nS
        % Assign common input fields, then override individually
        sIt = sIc;
        fn = fieldnames(sIuse(si));
        for fi = 1:length(fn)
            sIt.(fn{fi}) = sIuse(si).(fn{fi});
        end
        sI(si) = sIt;
    
        % Run scaling
         [dat(si).trackFits,dat(si).fitArrays,dat(si).in] = scaleTracks(tk,sI(si));
    end
    
    
    % Assign tracks and events for positioning and plotting (original order
    % was different from the manuscript events 1-2-3)
    whichEvents = {'25A4','24A','25B'};
    ne = length(whichEvents);

    % --------- NOTE --------
    % Tracks 2, 12, 20, 21 benefitted from additional constraints from mask track
    % widths to get good linear radius measures (good meaning over a longer height interval).
    % Those combined measures are contained in the sI2 run.
    % fits into main set (sI2 --> sI1) and to re-arrange sI9: uncomment
    % below
    datAll = dat;
%     dat(1).trackFits([tkNoOut]) = dat(2).trackFits;
    ff = fieldnames(dat(1).fitArrays); for fi=1:length(ff); dat(1).fitArrays.(ff{fi})(tkNoOut,:) = dat(2).fitArrays.(ff{fi}); end
    dat = dat(1);
    % -----------------------
    
    for si = 1:nS
        if save_scalingOutput
            save(fullfile(trackDir,sI(si).oFile),'dat','sI')
        end
    end
end

if flag_unsteady
    unsteadinessSandboxingV2
end


%% Params for plotting/saving scaling results figure

% Plot results
if plotResults
    [f1,f2]=plotScalingResultsMS(tk,dat(1),whichEvents);
end

if print_sampleTrackResults
%     figs = plotTrackFitsV2(tk(sI7.trackSet),dat(1),sI(1)); % Needs runs sI7
    [figs,ax] = plotTrackFitsV2(tk,dat(1),sI(1)); % Needs runs sI7
%     [figs2,ax2] = plotTrackFitsV2(tk,datAll(2),sI(2));
    for ai = 1:size(ax,1) 
        ll(ai)=letterlabel(['(' char(96+ai) ')'],ax(ai,1),12,[-0.23 1.01],[0 0 0],'bold'); %,'MyriadPro');
    end
    suppDims = [18 20];
    dpi = 500;
    for fi = 1:length(figs)
        figure(figs(fi))
        printpdf([suppTkFig '_' num2str(fi)],figDir,suppDims)
    end
end
if print_scalingResults
    oDims = [19 20];
    oDims2 = [9.5 11];
    [f1,f2]=plotScalingResultsMS(tk,dat(1),whichEvents);
    figure(f1)
    printpdf(scaleFig,figDir,oDims)
    figure(f2)
    printpdf(zBfig,figDir,oDims2)
end
