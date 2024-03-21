%% ++++++++++ EVENT 25A0 ++++++++++++++  
EVENT = '25A0'  ;  
run project_25A0

%         dataDir   = '~/Kahuna/data/sabancaya_5_2018/';
%         thermDir   = fullfile(dataDir,'image_exports/25A-0/');
%         matDir  = fullfile(thermDir,'reg-mat/');

        heads = matHeads; % fullfile(matDir,'frameHeads.mat');
        
        procDir = fullfile(thermDir,'preProc-mat/');

        % Use this with fram_tile to view plumeTrack output and iterate
        plumeTrackDir = fullfile(thermDir,'PTresults/');
        % pTfile  = fullfile(procDir,'PTresults/plumeTrack_output.mat'); % plumeTrack output file
        pTfile  = fullfile(plumeTrackDir,'plumeTrack_output.mat'); % plumeTrack output file
        % pTfile  = [];
        % T_scale   = 425;
        % dT_scale  = 400;
        % rel_scale = 0.4;

        % -------------------------------------------------------------------
        thresh_fg   = 276; % Temperature threshold applied to REFERENCE image to cut out foreground
        satVal      = 420; % Saturation temperature for thermal vid
        nullVal     = 190;    % Min value to scale temperatures down to

        T0          = 250;    % Middle value of smooth heaviside filter
        dT          = 10;     % Approximate 1/2-width of heaviside filter step
        kw          = 1/10;   % Width of filter...~2/dT?

        % Search for a histogram peak in these temperature ranges to visualize time
        % evolution. Range applied AT OR ABOVE given index
%         peakRange = [370 568 700 1576; ... % Index
%                      264 255 220 220; ... % Tmin
%                      310 310 310 310];    % Tmax
        peakRange = [];

        % This allows the thresholding values to be a changing function of index.
        % Leave empty to use a flat value
%         ITW = [378 420 550 600 700 800 1300 1567; ... % Index
%                267 267 257 254 249 246  240  240; ... % T0
%                 10  10  10  10  10  10   10   10];    % dT
        ITW = [];
        interp_meth = 'pchip';

        % thresh_fg = [217 240 260];
        bins     = nullVal:2:satVal;

        % -------------------INDICES-----------------------------
        ref_idx = 370;

        % Idx = 742;
        % Idx = [378 694:746]; %[378:4:500];
        % Idx = [378:520];
        % Idx = [ref_idx 550];% 450 700]; %SCRIPT ASSUMES FIRST INDEX IS THE REFERENCE IMAGE
        % Idx = [ref_idx 378:20:1500];
        % Idx = [383:10:600];
        % Idx = [543:50:1498];
        % Idx = [ref_idx 378:2:1576];
        Idx = [ref_idx 2:1000];
        % Idx = [ref_idx 553 853 1303 1503];
        % Idx = 1300;

        %%:400 402:5:450]; Idx = fliplr(Idx);

        %------ MANUAL POLYGONS to force plume mask or no mask? -----------
        design_polys = false;
        PolyFile     = fullfile(procDir,'manual_polygons_all.mat'); % Saves and loads polys from here
        % edit_polys   = 2; % If non-empty or non-zero, will attempt to load the Polygon with this index from PolyFile for editing?
        show_polys   = true; % Load above file and plot polies on frames

        npolys = [false]; % false]; % Boolean vector of length(number of polygons). true = force plume, false = force no plume;
        % IdxPoly = {[926:1448]}; % Indices for which to apply polygons (each cell is a list of indices for a 1 polygon)
        IdxPoly = {[1436:1576]};
        % Turn on to plot frames in DesPol_Idx and manually create polygons
        % DesPol_Idx = {[926 968 1268 1388]}; %~1-4 to plot and use for creating polies by click (1 row for each poly)
        DesPol_Idx = {[1436 1490 1538]};

        % ----------------PLOTTING FLAGS AND PARAMS--------------------
        design_flag = true; % Use this when designing the filters from a reference image - allows plotting in separate figures
        plot_flag   = true; % Use this when designing or vetting, turn off for mass production...or turn on and save?...
          imstat    = true; % Plot filter, histogram, and images steps for each file
          plot_idx  = [ref_idx 50 300 900 1488]; % Select indices to plot imstat\
          hist_tile = true; % Tile histograms for each file
          fram_tile = true; % Tile RAW frames...with outlines if provided
          plot_idx2 = [800:30:1500]; %[ 928:20:1576]; % Select indices for which to plot hist and frame tiles
        write_flag  = false; % Use this to write scaled and masked .mat frames

        % For hist tile plot and frame tile plot
        rows = 5;
        cols = 5;

        % Display position preferences to making image comparisons easier
        dd_filt = [50 800 600 400];
        dd_hist = [700 800 600 400];
        dd_mask = [1350 800 600 400];
        dd_rawi = [50 50 600 400];
        dd_mski = [700 50 600 400];
        dd_scli = [1350 50 600 400];