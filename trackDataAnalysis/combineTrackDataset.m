% THIS SCRIPT IS NOW OUTDATED AND REPLACED with buildTrackDataset in folder
% pulseTrackDataAnalysis/

% Build final, clean data sets
clear all; close all
%% Input data
dataDir     = '~/Kahuna/data/sabancaya_5_2018/image_exports'; %';

% Excel file with data selection parameters
analysisXLSX = '~/Kahuna/data/sabancaya_5_2018/Data_and_Event_index.xlsx';
analysisSheet = 'pulseTrack_notes';

oFile = sprintf('~/Kahuna/data/sabancaya_5_2018/pulseTrack_analysis/allTracks_raw_24A_25A4_25B_%s.mat',datestr(now,'YYYY-mm-dd'));
% ---- 24A -----
iD24A = 1;
iDat(iD24A).event = '24A'; % Event name
iDat(iD24A).thermCube = fullfile(dataDir,'/24A/thermCubeAnalysis/thermStats_2020-06-29_z644_x578_t1195_revCorr.mat');
iDat(iD24A).velCube   = fullfile(dataDir,'/24A/thermCubeAnalysis/opticFlowCNL_20-06-29_n1195_filtered.mat');
iDat(iD24A).trackFile = fullfile(dataDir,'/24A/thermCubeAnalysis/Tracks_24A_processed.mat');

iDat(iD24A).satVal = 424.74; 
% 
% plotThermVel
iDat(iD24A).Trange = [0 80];
iDat(iD24A).proi = [98 1400 -430 440];



% ---- 25A -----
iD25A = 2;
iDat(iD25A).event = '25A4'; % Event name
iDat(iD25A).thermCube = fullfile(dataDir,'/25A-4/thermCubeAnalysis/thermStats_2020-07-20_z704_x624_t2439.mat');
iDat(iD25A).velCube   = fullfile(dataDir,'/25A-4/thermCubeAnalysis/opticFlowCNL_20-07-31_n2439_filtered.mat');
iDat(iD25A).trackFile = fullfile(dataDir,'/25A-4/thermCubeAnalysis/Tracks_25A4_processed.mat');

iDat(iD25A).satVal = 420; % Not accurately known? 
% 
% plotThermVel
% cmap = gray(200);
iDat(iD25A).Trange = [0 35];
iDat(iD25A).proi = [98 1050 -500 150];
% iDat(iD25A).tracksI = 1:9;

% ---- 25B -----
iD25B = 3;
iDat(iD25B).event = '25B'; % Event name
iDat(iD25B).thermCube = fullfile(dataDir,'/25B/thermCubeAnalysis/thermStats_2020-10-26_z710_x571_t1289.mat');
iDat(iD25B).velCube   = fullfile(dataDir,'/25B/thermCubeAnalysis/opticFlowCNL_20-06-18_n1289_filtered.mat');
iDat(iD25B).trackFile = fullfile(dataDir,'/25B/thermCubeAnalysis/v2Tracks/Vtracks_25B_processed.mat');

iDat(iD25B).satVal = 403.1;

% plotThermVel
iDat(iD25B).Trange = [0 140]; % Plot T range
iDat(iD25B).proi = [98 1400 -430 440]; % Region of interest

% Start and end indices to use for each track


% ---- thingumits ----
% plot params
pdx = 0.05;
pdy = 0.09;
ppads = [0.07 0.03 0.1 0.03];

cmap = gray(200); % Plot colormap

figdir = '~/Kahuna/data/sabancaya_5_2018/image_exports/25B/figures';

% ----- TRACK STATS RETRIEVAL PARAMS -----


%% Run through all input Data, output combined, labeled tracking data set\

allTracks = [];
allTrkPar = [];

tcnt = 0;
for kk = 1:length(iDat)
    fprintf('\nEvent: %s\n',iDat(kk).event)

    load(iDat(kk).trackFile)

    allTrkPar = [allTrkPar trackPar];
    
    nTracks = length(Vtrack);
    for ii = 1:nTracks % Within event track count
        tcnt = tcnt+1; % Total track count
        
        % Copy initial fields
        dat(tcnt).event      = iDat(kk).event;
        dat(tcnt).eventTrack = ii;
        dat(tcnt).thermCube  = iDat(kk).thermCube;
        dat(tcnt).velCube    = iDat(kk).velCube;
        dat(tcnt).trackFile  = iDat(kk).trackFile;
        
        dat(tcnt).satVal     = iDat(kk).satVal;
%         dat(tcnt).   = Vtrack(ii).;
        dat(tcnt).imSz       = Vtrack(ii).cubeDims(1:2);
        dat(tcnt).zI         = Vtrack(ii).zI;
        dat(tcnt).tI         = Vtrack(ii).tI;
        dat(tcnt).t          = Vtrack(ii).t;
        dat(tcnt).clustZ     = Vtrack(ii).clustZ;
        dat(tcnt).clustIdx   = Vtrack(ii).clustIdx;
        dat(tcnt).clustI     = Vtrack(ii).clustI;
        dat(tcnt).npx        = Vtrack(ii).npx;

    end
    
end

tk = dat; clear dat;
fprintf('Saving:\n\t%s\n',oFile)
save(oFile,'tk','allTrkPar')
