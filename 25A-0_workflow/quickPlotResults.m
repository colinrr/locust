%% =====================================================================
%       25A-0 Quick results plot for thermal source history
% ======================================================================

clear ; close all;
% Set root data directories here
dataDir = fullfile('C:\Users\crowell\','/Kahuna/data/sabancaya_5_2018/image_exports/25A-0/thermCubeAnalysis');

% Core data files:
% Thermal source times series data
thermalSourceFile   = fullfile(dataDir,'thermSource_25A0_2024-03-20.mat'); 
% Imagery data cube structure
thermCube           = fullfile(dataDir,'thermStats_2024-03-19_z709_x618_t1136.mat'); 

%% Load data 
% Load the core data files
load(thermCube)
load(thermalSourceFile)

%% Run plots - SEE ALSO plots at the bottom of 'workflowDriver_25A0.m'
plotFrames = 255; % Targeted frame for second figure

% Plot params
dx = 0.1;
dy = 0.13;
ysz = [0.5 1];
% Plot the source time series - no script dependencies
figure('position',[100 200 1000 900])
ax(1) = tightSubplot(2,1,1,dx,dy,[],[],ysz);
plot(S.t_since_00utc,S.dT95, S.t_since_00utc, S.dTmax,'LineWidth',1.3)
  % '-> see S.t_local, S.t_utc, S.t for equivalent time series
  
hold on
xlabel('seconds since 00:00 UTC')
ylabel('Excess Temp., \Delta T (K)')
grid on
axis tight
plot(S.t_since_00utc(plotFrames)*[1 1], ylim,'--k','LineWidth',1.5)
legend('95th percentile','Maximum','Plotted frame')

% Plot a closeup of the source region for a few frames
zoomRegion = [D.z(1) 300 -400 -50];
deltaTempRange  = [-30 40];

ax(2) = tightSubplot(2,2,3,dx,dy,[],[],ysz);  % Full frame
% figure('position',[100 100 1200 600],'Name','Source retrieval region')
plotThermVelocities(D.x,D.z,[],[],'Mask',D.sourceMask,'Thermal',D.T,...
    'Trange',deltaTempRange,'idx',plotFrames,...
    'atmo',D.atmo,'time',S.t_since_00utc)
axis tight

ax(3) = tightSubplot(2,2,4,dx,dy,[],[],ysz); % Source zoom
plotThermVelocities(D.x,D.z,[],[],'Mask',D.sourceMask,'Thermal',D.T,...
    'Trange',deltaTempRange,'idx',plotFrames,...
    'atmo',D.atmo,'ROI',zoomRegion,'time',S.t_since_00utc)
