% Keeping Irbis processing separate from pure matlab image analysis
%
% ISSUE: irb2ascii program, when applied to and Irbis VIDEO, outputs
%  correct meta-data, but just a matrix of zeros for the body. Conversely,
%  exporting the frames individually from IRT_analyzer gives the full body,
%  but incomplete headers (specifically the timestamps are not recorded to
%  the millisecond, making any real-time processing impossible). On top of
%  which, IRT_analyzer output will drop frames, so we have to figure out
%  which frames were dropped without the benefit of accurate time stamps.
%  Yay.
%
% SOLUTION: (annoyingly) 
%  1) Export images in irbis2ascii twice - once from raw
%  video file, once from IRT_analyzer. Read meta-data from raw export, body
%  info from IRT analyzer. SO! Get the corresponding folders correct!
%    IRT_analyzer OUTPUT --> BODY (aka Image) data
%    IRB direct OUTPUT   --> HEADER data
%
% * Where able to export from IRBIS program orginally, this is not an issue
% and can go direct to irbAsc2Mat.
%
clear all; close all
% Set path
% homeDir = '~';
% addpath(genpath(fullfile(homeDir,'code/research-projects/sabancaya/pulseTracker-dev/IRBtools')))

homeDir = 'C:\Users\crowell\';
addpath(genpath(fullfile(homeDir,'\Documents\GitHub\pulseTracker-dev\')))
%%

% ASCII directory
ascDir  = fullfile(homeDir,'Kahuna\data\sabancaya_5_2018\25_may_2018_morning\JENOPTIK_IR_PICKLES\180525AA\');
% ascDir  = fullfile(homeDir,'\Kahuna\data\sabancaya_5_2018\25_may_2018_afternoon\180525BI\');

dataDir   = fullfile(homeDir,'/Kahuna/data/sabancaya_5_2018/');
% -------- Set up directories and files for this event --------
% 25A-0
matDir   = fullfile(dataDir,'image_exports/25A-0/mat/');
ascDir_body = fullfile(ascDir,'25A-0-body-ascii/');
ascDir_head = fullfile(ascDir,'25A-0-head-ascii/');


IRBheads = fullfile(matDir,'frameHeads_crop_to_dropped_frames.mat'); % Good timestamps - best guess was that 3 frames were dropped so N=1139
IRTheads = fullfile(matDir,'frameHeads_IRT_bad_timestamps.mat'); % Bad timestamps attached to data, N=1136

% CORRECTED output headers
frameHeads = fullfile(matDir,'frameHeadsFixed.mat');

% glob_spec = '*.tif';
% glob_spec = 'RT_1030_corrected_*.txt';
glob_spec = '*.txt';

% paramsIRB = fullfile(matDir,'params-heads.mat');

%% ================= DRIVER SWITCHES AND WORKFLOW ==================
%  ------ Data conversion and from IRBIS ascii -------
% flag_tif2mat    = false; % Out of date, use asc2mat
flag_readHeads   = false;
flag_asc2mat     = false;  % IRB Ascii to .mat conversion
% addTLframes25B %--> DO NOT RUN EXCEPT ON A FRESH OUTPUT OF asc2mat FILES
flag_fixTimes    = true;

%% ========================== DO THE THING ==========================

if flag_readHeads
    [Th,IRBheads] = irbHeads2Table(ascDir_head,matDir,glob_spec);
end
if flag_asc2mat
%     [Tb,IRTheads] = irbAsc2Mat(ascDir_body,matDir,IRBheads,glob_spec);
    [Tb,IRTheads] = irbAsc2Mat(ascDir_body,matDir,'None',glob_spec);
end
if flag_fixTimes
    [Tfix,Icut] = fixTableTimes(IRBheads, IRTheads, matDir);
end

%% Quick check for corrected timestamps
if flag_fixTimes
    T_head_all = load(IRBheads);
    T_bad = load(IRTheads);
    T_fixed = load(frameHeads);
    
    figure
    % show time vectors
    subplot(2,1,1)
    plot(T_bad.T.Time,'Color',[0.7 0.7 0.7]);
    hold on
    plot(T_head_all.T.Time);
    plot(T_fixed.T.Time);
    
    % Show time diff
    subplot(2,1,2)
    plot(diff(T_bad.T.Time),'Color',[0.7 0.7 0.7]);
    hold on
    plot(diff(T_head_all.T.Time));
    plot(diff(T_fixed.T.Time));
end