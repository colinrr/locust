%% ================= PARAMETER INPUT ==================

clearvars -except D; % close all;
%---- Build path ----
addpath(genpath('../.'))

% ==== SELECT PROJECT =====
run project_25A0 % Load project directories & parameters

if ~exist(homeDir,'dir')
    error('Check project directories!')
end

%% ================= DRIVER SWITCHES AND WORKFLOW ==================
% √ = DONE; X = skipped
%  ------ Data conversion and registration -------
flag_pixelreg    = false;  % X  Register thermal images (correct shaking etc)

% ------ plumeTracker (mask generation), geometry and initial calcs ------
flag_mapPixels   = false;   % √ Generate mapping function to convert pixels to meters
flag_preProcess  = true;   % Copy and pre-process thermal images for plumeTracker
    get_foreground_mask   = false; % √ Design a single mask file to remove foreground
    view_histograms       = false; % √ Useful for designing filters
    % test_filters   = 
    apply_filters_test    = false; % √  Apply "filters" to test frames
    artificial_ref        = false; % X  Artificially zero out reference image to nullVal
    flag_plumetrack_test  = false; % √ Test first 30 frames with plumeTracker
    apply_filters         = false; % √ Apply filters to all frames
flag_plumetrack  = false;   % √ Run Bombrun plume segmentation algorithm

flag_poly2mask   = false;   % X Apply manual polygons to plume masks
flag_plumecalcs  = false;   % X Basic H, v, A calcs for segmented plumes

% ----- Data Cube Workflow: re-grid images, get thermCube, atmoProfile --------
flag_interpTherm = false;   % √ Interpolate thermal images/masks to regular x,z grids
flag_thermCube   = false;   % √ Get thermal data cube
flag_atmoProfile = false;    % √ Get atmospheric temperature profiles

% ----- Data Cube Workflow: velocimetry --------
flag_opticFlow    = false;   % Optical Flow analysis on thermal or gradient video
flag_filtVelocity = false;   % Apply low-pass filter to optic flow velocities
flag_timeAverage  = false;   % Get time-averaged plume image and profiles

% ----- Data Cube Workflow: pulseTracker --------
flag_getSource    = false;   % √ Run source pulse detection, get source history
flag_pulseTrack   = false;  % Track detected sources.
flag_postProcess  = false;  % Apply smoothing and pixel exclusion to tracks

% ---- Analysis ----


% ------ pretty pictures --------
flag_image2dem   = false;   % Project thermal images onto the DEM
flag_riseDiag    = false;   % Plot rise diagram

% ------ Deprecated Flags -----
% flag_spectra1D
% flag_gradientVid  = false;   % Make video of thermal gradient from thermal data cube
% flag_scaledVid    = false;   % Make scaled thermal video for optic flow analysis
% flag_thermCorr    = false;   % Cross-correlation analysis on thermal data cube

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
    preProcWorkflow_25A0
end
if flag_plumetrack
    %                                                             oDir,      ref, deb, fin, dt
    [content,refData] = mainTrackPlume(ptInputDir,ptInputHeads,plumeTrackDir, ref, deb,fin,dN,fixpix,PTplotFlags,delT); 
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
    interpThermal(interpIn,interpDir,{matHeads,ptHeads},geomf,interpIdx,[],vmax,polyFile)
    
    % Re-grid reference images for atmospheric profile retreival
%     interpThermal(interpIn,interpRefDir,{regHeads,ptHeads},geomf,interpRefIdx,[],vmax,polyFile)
    % Long, decimated time series
%     interpThermal(interpIn,interpRefDir2,{regHeads,ptHeads},geomf,interpRefIdx2,[],vmax,polyFile)
end
if flag_thermCube
%     [D,thermCube] = getThermCube(interpDir,interpHeads,geomf,thermIdx,ROI,cubeDir);
    [D,atmoRefCube] = getThermCube(interpDir,interpHeads,geomf,thermIdxAtmo,ROI,cubeDir); % Reference cube for atmo profile
end
if flag_atmoProfile
    % PART 1 - retrieve satellite profiles and get initial fit
%      [atmo,atmoFile] = getAtmoProfile(atmoHDF,vent(1:2),localTimezoneOffset); % Viewing atmo profiles only

      % Comparing w/ thermal data
%     [atmo,initialStats,atmoFile] = getAtmoProfile(atmoHDF,vent(1:2),localTimezoneOffset,atmoRefCube,atmoRefFrames,atmoTthresh); % Comparing w/ thermal data, no save
    [atmo,initialStats,atmoFile] = getAtmoProfile(atmoHDF,vent(1:2),localTimezoneOffset,atmoRefCube,atmoRefFrames,atmoTthresh,thermDir); % Comparing w/ thermal data, save output

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
%     loadif(velCube,'V')
    % Avoid the need for optic flow for now 
    loadif(thermCube,'D')
%     S = getThermSource(D,[],thermSource); %trackParams,detParams,detection_plotflag);
    
    % Save output
    S = getThermSource(D,[],thermSource,'outputFile',thermSourceFile);
    
    % Adding some extra tidbits to the data files - datetime vectors &
    % source region mask
    % Pull out T95, Tmax, t into column vectors for ease of use
    load(interpHeads)
    S.dTmax = S.T0.max';
    S.dT95  = S.T0.prctile(5,:)';
    S.t = S.T0.t';
    S.t_local = datetime(T.Timestamp(1),'TimeZone','America/Lima',...
        'Format','y-MMM-d HH:mm:ss.SSS Z',"ConvertFrom","datenum") + seconds(S.t);
    S.t_utc   = datetime(S.t_local,'TimeZone','UTC');
    utc_00 = datevec(S.t_utc(1)); utc_00 = datetime([utc_00(1:3) 0 0 0],'Format','y-MMM-d HH:mm:ss.SSS Z','TimeZone','UTC');
    S.t_since_00utc = seconds(S.t_utc - utc_00);

    save(thermSourceFile,'S','-append')
    for ii=1:size(D.mask,3) 
        [~,D.sourceMask(:,:,ii),D.sourcePoly(ii)] = getROI(D.mask(:,:,ii),'iLims',[1 20],'maxRegions',1); 
    end
    save(thermCube,'D','-append')

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
% playFrames(matDir,matHeads,15,[1:10:1136],[],[],[230 330]) % This will play frames ~1 s apart through whole sequence

% regFrameDiff(matDir,matHeads,regDir,regHeads,regROI,refIdx)
% [700:790 1020:1140]

% >> Calculate mean or max difference between pixels in an roi for all frames
% >> in 'idx'. Useful to show which frames may need registration.
% load(regParams)

% Used these routines -vv- to decide registration was not needed:
% [Fmean,Fmax,didx]=regFrameDiff(matDir,matHeads,[],[],regROI); % Pre-registered
% playFrames(matDir,matHeads,5,[563:618],[],regROI([3 4 1 2]),[230 330]) % This will play frames ~1 s apart through whole sequence
% playFrames(matDir,matHeads,5,[670:690],[],regROI([3 4 1 2]),[230 330]) % This will play frames ~1 s apart through whole sequence


%------- plumeTracker (mask generation), geometry and initial calcs -------
 % Play single image to check reference pixel location
% playFrames(matDir,matHeads,[],1,[],[],[230 330])

% plotImageMap(geom,fullfile(matDir,'\25A_00_0001.mat')) % Plot image projection

% playFrames(interpDir,interpHeads,15,interpIdx,[],[],[230 350],true)
% '-> ALSO VIEW VIDEO/IMAGE OUTPUT in plumeTracker dir (usually plumeTrackDir)

% ------ Check interpolated or plumeTrack frames with masks ------
% playFrames(interpDir,interpHeads,6,interpIdxTest,[],[],[230 350],true)
% vv Compare a frame projection pre- and post- re-gridding in interpTherm vv
% frameProjectionPair({matDir;matHeads},geomf,{interpDir;interpHeads},[],ref_idx,'stack')

% ------ Check thermCube output (getThermCube) -------
% loadif(thermCube,'D');
% figure
% plotThermVelocities(D.x,D.z,[],[],'Mask',D.mask,'idx',round(linspace(1,length(D.t),100)),'Thermal',D.T)

% ------ Check atmo profile fit (getAtmoProfile) -------
% '-> show example of cropped plume mask used for fitting
% [~,cropMask,~] = getROI(D.mask(:,:,1),'iLims',atmoImgROI(1:2),'jLims',atmoImgROI(3:4));
% plotThermVelocities(D.x,D.z,[],[],'Mask',cropMask,'idx',1,'Thermal',D.T)

% '-> show the video sequence with removed atmo profile and tightened mask
% plotThermVelocities(D.x,D.z,[],[],'Mask',D.mask,'idx',1:10:1136,'Thermal',D.T,'atmo',D.atmo,'Trange',[-30 30])

% ------ Check optical flow output (thermOpticFlow) -------
% plotThermVelocities(V.x,V.z,V.Vx,V.Vz,'Mask',D.mask(:,:,1:end-1),'idx',[],'Vmax',20)
% plotThermVelocities(V.x,V.z,V.Vx,V.Vz,'Mask',D.mask,'idx',1:20:1100,'Vmax',20,'Thermal',D.T)

% ------ Check thermal source history output (getThermSource) ------
% '- > Plot temperature statistics for a retrieved source history (from
% getThermSource)
[ax,ll]=plotThermStats(S.T0,[],'label','\Delta T','unit','K'); % Temperature

% '-> Show video sequence highlighting source region retrieval
figure('position',[100 100 1000 800],'Name','Source region retrieval')
plotThermVelocities(D.x,D.z,[],[],'Mask',D.sourceMask,'Thermal',D.T,...
    'Trange',[-30 30],'idx',1:5:1136,...
    'atmo',D.atmo,'ROI',[D.z(1) 750 D.x(1) 500],'time',D.t)

