%% ================= PROJECT INPUT ==================


%---- Build path ----
% addpath(genpath('thermImagePreprocessing'))
% addpath(genpath('plumeTracking'))
% addpath(genpath('sabancayaScripts'))
% addpath(genpath('featureTracking'))
% addpath(genpath('utils'))
addpath(genpath('.'))
% addpath(genpath('ijcv_flow_code')) % For optical flow toolbox

% -------------------

% ==== SELECT PROJECT =====
% - These project files set all working directories and input parameters
% for each event -

% run sabancayaScripts/project_event1
% run sabancayaScripts/project_event2
run sabancayaScripts/project_event3


if ~exist(homeDir,'dir')
    error('Check project directories!')
end

%% ================= DRIVER SWITCHES AND WORKFLOW ==================
% FLIP the switches below to run different elements of the workflow, and
% don't forget to check the plotting routines down at the bottom.

%  ------ Data conversion and registration -------
flag_pixelreg    = false;  % Register thermal images (correct shaking etc)

% ------ plumeTracker (mask generation), geometry and initial calcs ------
flag_mapPixels   = false;   % Generate mapping function to convert pixels to meters
flag_preProcess  = false;   % Copy and pre-process thermal images for plumeTracker
flag_plumetrack  = false;   % Run Bombrun plume segmentation algorithm
flag_poly2mask   = false;   % Apply manual polygons to plume masks
flag_plumecalcs  = false;   % Basic H, v, A calcs for segmented plumes

% ----- Data Cube Workflow: re-grid images, get thermCube, atmoProfile --------
flag_interpTherm = false;   % Interpolate thermal images/masks to regular x,z grids
flag_thermCube   = false;   % Get thermal data cube
flag_atmoProfile = false;    % Get atmospheric temperature profiles

% ----- Data Cube Workflow: velocimetry --------
flag_opticFlow    = false;   % Optical Flow analysis on thermal or gradient video
flag_filtVelocity = false;   % Apply low-pass filter to optic flow velocities
flag_timeAverage  = false;   % Get time-averaged plume image and profiles

% ----- Data Cube Workflow: pulseTracker --------
flag_getSource    = false;   % Run source pulse detection, get source history
flag_pulseTrack   = false;  % Track detected sources.
flag_postProcess  = false;  % Apply smoothing and pixel exclusion to tracks

% ---- Analysis ----


% ------ pretty pictures --------
flag_image2dem   = true;   % Project thermal images onto the DEM
flag_riseDiag    = false;   % Plot rise diagram

%% ========================== DO THE THING ==========================
%  ------ Image registration/stabilization -------
if flag_pixelreg
    [regParams,regHeads] = thermPixelReg(matDir,matHeads,regTscale,regIdx,refIdx,regROI,regDir);
end

% ------ plumeTracker (mask generation), geometry and initial calcs ------
if flag_mapPixels
    [geom,geomf] = mapPixels(obsLLE,vent,refLLE,refPix,hfov,vfov,imsz,regParams,thermDir);
end
if flag_preProcess
    preProcThermal
end
if flag_plumetrack
    %                                                             oDir,      ref, deb, fin, dt
    [content,refData] = mainTrackPlume(ptInputDir,ptInputHeads,plumeTrackDir, ref, deb,fin,dN,fixpix,PTplotflags,delT); 
end
if flag_poly2mask
    [T,Ref,update_time] = maskPolyClip(ptHeads,polyFile,true);
end
if flag_plumecalcs
%     load(geomf)
    [dat,geom] = plumeTrackCalcs(ptHeads,obsLLE,refLLE,refPix,geomf,plotCalcs);
end

% ----- Data Cube Workflow: re-grid images, get thermCube, atmoProfile --------
if flag_interpTherm
    % Re-grid main image sequence
    interpThermal(interpIn,interpDir,{regHeads,ptHeads},geomf,interpIdx,[],vmax,polyFile)
    
    % Re-grid reference images for atmospheric profile retreival
%     interpThermal(interpIn,interpRefDir,{regHeads,ptHeads},geomf,interpRefIdx,[],vmax,polyFile)
    % Long, decimated time series
%     interpThermal(interpIn,interpRefDir2,{regHeads,ptHeads},geomf,interpRefIdx2,[],vmax,polyFile)
end
if flag_thermCube
    [D,thermCube] = getThermCube(interpDir,interpHeads,geomf,thermIdx,ROI,cubeDir);
end
if flag_atmoProfile
    % PART 1 - retrieve satellite profiles and get initial fit
%      [atmo,atmoFile] = getAtmoProfile(atmoHDF,vent(1:2),localTimezoneOffset); % Viewing atmo profiles only

      % Comparing w/ thermal data
%     [atmo,initialStats,atmoFile] = getAtmoProfile(atmoHDF,vent(1:2),localTimezoneOffset,atmoRefCube,atmoRefFrames,atmoTthresh); % Comparing w/ thermal data, no save
%     [atmo,initialStats,atmoFile] = getAtmoProfile(atmoHDF,vent(1:2),localTimezoneOffset,atmoRefCube,atmoRefFrames,atmoTthresh,thermDir); % Comparing w/ thermal data, save output

    % PART 2 - run final profile fit on filtered data
%     [D,atmo,fitStats] = fitAtmoProfile(atmoRefCube,atmoFile,atmoRefFrames,atmoImgROI,TerrorThresh,zUncertaintyThresh,atmoTthresh); % Atmo from file (writes updated file) 
    [D,atmo,fitStats] = fitAtmoProfile(atmoRefCube,atmo,atmoRefFrames,atmoImgROI,TerrorThresh,zUncertaintyThresh,atmoTthresh); % Atmo from workspace(no save to disk)
end


% ----- Data Cube Workflow: velocimetry --------
if flag_opticFlow
    [V,velCube] = thermOpticFlow2O(thermCube, opticParams, 'opticFlowCNL');
end
if flag_filtVelocity
    loadif(thermCube,'D')
    [V,Vspec] = filterVelocities(velCube,D.mask,[],OFfiltlim,[],filtTest,filtPlot); % write mode
%     [V,Vspec] = filterVelocities(V,D.mask,[],OFfiltlim,[],false,true); % non-write mode
end
if flag_timeAverage
    loadif(velCube,'V')
    loadif(thermCube,'D')
%     [Tavg,maskAvg] = getAveragedImage(D.T,D.mask,[],[],avgIdx,D.atmo); % No velocities
    Tavg = getAveragedImage(D.T,D.mask,V.Vx,V.Vz,D.x,D.z,avgROI,avgIdx,D.atmo,avgWinHeight,D.geom); % with Vx,Vz
    % Manual save...
    save(avgImg,'Tavg')
end


% ----- Data Cube Workflow: pulseTracker --------
if flag_getSource
    loadif(velCube,'V')
    loadif(thermCube,'D')
    S = getThermSource(D,V,thermSource); %trackParams,detParams,detection_plotflag);
%     S = getThermSource(D,V,thermSource,'outputFile',thermSourceFile); % Save output
end
if flag_pulseTrack
    loadif(velCube,'V')
    loadif(thermCube,'D')
    
    % -----
    % Structure tracking can be run directly from here using the project_event
    % scripts, but in practice was run separetely using the
    % "sabancayaScripts/structureTracking" scripts below.
%     [Vpos,trackPar] = trackStructure(V,D,varargin{:}); 

%     run sabancayaScripts/structure_Tracking_event1
%     run sabancayaScripts/structure_Tracking_event2
    run sabancayaScripts/structureTracking_event3
    % -----
end
if flag_postProcess
    loadif(thermCube,'D')
    load(rawTrackFile)
    Vtrack = postProcessTracks(Vtrack,trackPars,D.maskSmoothLength,trackFile);
end



%% =======================================================================
%                   PulseTracker PLOTTING TOOLS
% =======================================================================
% Use the functions below for workflow QC and visualization!

if flag_image2dem
    dem_frames = 2; %[1:size(dat,1)];
    lp=image2dem(dat,px2m,geom,demfile,dem_roi,matDir,dem_frames);
end
if flag_riseDiag
%     D = plotRiseDiagram(matDir,fullfile(outputDir,'geometry.mat'),ventPix,profileType,Ridx);
end

%--------------- Data conversion and registration ---------------- 
% > Playthrough of .mat frames for quick viewing
% playFrames(regDir,regHeads,15,interpIdx(1:5:end),[],[],[230 400])  % Registered frames
% playFrames(interpDir,interpHeads,15,interpIdx(1:5:end),[],[],[230 400]) % Registered and re-gridded frames

% > Calculate mean or max difference between pixels in an roi for all frames
% > in 'idx'. Useful to show which frames may need registration.
% load(regParams)
% [Fmean,Fmax,didx]=regFrameDiff(matDir,matHeads,[]); % Pre-registered
% [Fmean,Fmax,didx]=regFrameDiff(matDir,matHeads,regDir,regHeads,R.newROI); % Compare pre- and post-registered

%------- plumeTracker (mask generation), geometry and initial calcs -------
% -- Plot image projection---
% load(geomf)
% plotImageMap(geom) 

% -- Plot image projection with sample image---
% load(geomf)
% load(regHeads) % Get ref image for plotImageMap
% plotImageMap(geom, fullfile(regDir,T.File{num2str(interpIdx(1))}))
% '-> ALSO VIEW VIDEO/IMAGE OUTPUT in plumeTracker dir (usually plumeTrackDir)

% ------ Check interpolated or plumeTrack frames with masks ------
% vv Compare a frame projection pre- and post- re-gridding in interpTherm vv
% frameProjectionPair({regDir;regHeads},geomf,{interpDir;interpHeads},[],interpIdx(1),'stack')

% ------ Check thermCube output (getThermCube) -------
% loadif(thermCube,'D');
% loadif(velCube,'V');
% figure
% plotThermVelocities(D.x,D.z,[],[],'Mask',D.mask,'idx',round(linspace(1,length(D.t),50)),'Thermal',D.T)

% ------ Check optical flow output (thermOpticFlow) -------
% loadif(thermCube,'D');
% loadif(velCube,'V');
% plotThermVelocities(V.x,V.z,V.Vx,V.Vz,'Mask',D.mask,'idx',[],'Vmax',20,'roi',[100 1500 -500 600])  % Velocity only
% plotThermVelocities(V.x,V.z,V.Vx,V.Vz,'Mask',D.mask,'idx',[],'Vmax',20,'Thermal',D.T,'roi',[100 1500 -500 600])  % Velocity and thermal



% ------ Check thermal source history output (getThermSource) ------
% -- > Plot temperature statistics for a retrieved source history (from
% getThermSource)
% ax=plotThermStats(S.Td,[],S.tTrig);
% title(ax(1),'Source Statistics, Detection Window')
% ax=plotThermStats(S.T0,[],'tTrig',S.tTrig); % Temperature
% ax=plotThermStats(S.V0,[],'tTrig',S.tTrig); % Velocity
% title(ax(1),'Source Statistics, Tracking Window')
 % --> Show detection windows against a Rise Diagram
% plotRiseDiagram(D.z,D.t,D.T,'detections',S) % Basic
% plotRiseDiagram(D.z,D.t,D.T,'mask',D.mask,'atmo',D.atmo,'detections',S) 

