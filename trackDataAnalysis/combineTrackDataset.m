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
% cmap = gray(200);
iDat(iD24A).Trange = [0 80];
iDat(iD24A).proi = [98 1400 -430 440];
% iDat(iD24A).tracksI = 1:9;

% iDat(iD24A).tcutT = []
%     1 370 % 1
%     22 70 % 2
%     6 126 % 3
%     1 116 % 4a
%     1 123 % 4b
%     1 177 % 5a
%     1 168
%     1 111
%     1 158
%     ];

% iDat(iD24A).tcutR = [];

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

% iDat(iD25B).idx = [150 200 285 320]; % 350]; %300; % 285
iDat(iD25B).satVal = 403.1;

% plotThermVel
iDat(iD25B).Trange = [0 140]; % Plot T range
iDat(iD25B).proi = [98 1400 -430 440]; % Region of interest
% iDat(iD25B).tracksI = [1 2 3 4 5];  % Which tracks to keep

% Start and end indices to use for each track
% iDat(iD25B).tcutT = [
%     33 425  % 1
%     84 200  % 2
%     117 305 % 3
%     70 428 % 4
%     97 640]; % 5

% iDat(iD25B).tcutR = [];

% ---- thingumits ----
% plot params
pdx = 0.05;
pdy = 0.09;
ppads = [0.07 0.03 0.1 0.03];

cmap = gray(200); % Plot colormap

figdir = '~/Kahuna/data/sabancaya_5_2018/image_exports/25B/figures';

% ----- TRACK STATS RETRIEVAL PARAMS -----


%% Run through all input Data, output combined, labeled tracking data set\
% co0 = get(gca,'ColorOrder');
% co0 = co0([1:2 4:7],:);
% co = repmat(co,[ceil(nTracks/size(co,1)) 1]);
% co = co(1:nTracks,:);
% dat = iDat;
allTracks = [];
allTrkPar = [];

tcnt = 0;
for kk = 1:length(iDat)
    fprintf('\nEvent: %s\n',iDat(kk).event)
%     clear D V
%     loadif(iDat(kk).thermCube,'D')
%     loadif(iDat(kk).velCube,'V')
    load(iDat(kk).trackFile)
%     disp('Clipping mask...')
%     [~,maskClip] = removeAtmoProfile(D.mask,D.T,D.atmo,trackPar(1).nHalfMax,false);

    allTrkPar = [allTrkPar trackPar];
    
%     nTracks = length(iDat(kk).tracksI);
%     textprogressbar('Track Stats:  ')
%     textprogressbar(0)
    nTracks = length(Vtrack);
    for ii = 1:nTracks % Within event track count
%         ti = iDat(kk).tracksI(ii);
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
%         dat(tcnt).tcutT  = iDat(kk).tcut(ii);

        
%         dat(tcnt).Tstats = getThermStats(D.T,D.mask,D.z,D.t,Vtrack(ii).tI,Vtrack(ii).clustIdx,[],D.atmo,iDat(kk).satVal);
        % Absolute value for U?
%         dat(tcnt).Ustats = getThermStats(V.Vx,maskClip,D.z,D.t,Vtrack(ii).tI,Vtrack(ii).clustIdx,[]);
%         dat(tcnt).Wstats = getThermStats(V.Vz,maskClip,D.z,D.t,Vtrack(ii).tI,Vtrack(ii).clustIdx,[]);
        
%         dat(tcnt).t = Vtrack(ti).t(dat(tcnt).tcut(ii,1):dat(tcnt).tcut(ii,2)) - Vtrack(ti).t(dat(tcnt).tcut(ii,1));
        
        % Want a flag for plumeTracking mode....?
%         dat(tcnt).diam = Vtrack(ti).npx(tcut2(ii,1):tcut2(ii,2)).^.5*4/pi*D.dz;
%         dat(tcnt).diam0 = Vtrack(ti).npx(1).^.5*4/pi*D.dz;    
%         textprogressbar(ii/nTracks*100)
    end
%     textprogressbar(' --> Done')
    
end

tk = dat; clear dat;
fprintf('Saving:\n\t%s\n',oFile)
save(oFile,'tk','allTrkPar')
