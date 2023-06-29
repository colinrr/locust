% Compile T-averaged event images
clear all; close all

dDir = '~/Kahuna/data/sabancaya_5_2018/image_exports/';
avgImgFiles = {
               fullfile(dDir,'25A-4/thermCubeAnalysis/tAvg_25A4_steady_2022-07-16.mat')
               fullfile(dDir,'24A/thermCubeAnalysis/tAvg_24A_quasiSteady_2022-07-16.mat')
               fullfile(dDir,'25B/thermCubeAnalysis/tAvg_25B_allLargePulses_2022-07-16.mat')
               fullfile(dDir,'25B/thermCubeAnalysis/tAvg_25B_all_t0-87_2022-07-16.mat')
               fullfile(dDir,'25B/thermCubeAnalysis/tAvg_25B_early_t0-27_2022-11-02.mat')
               fullfile(dDir,'25B/thermCubeAnalysis/tAvg_25B_decay_t27-85_2022-07-16.mat')
               };
            
events = {'Event 1','Event 2','Event 3','Event 3','Event 3','Event 3'};
desc   = {'25A-4 Steady','24A Quasi-Steady Period','25B Large Pulses','25B Full History','25B Main Pulse','25B Decay Period'};
% added on because I was too dumb to do it before
timeSpans = {[0 308.25],[47.08 150.00],[0 54.33],[0 84.73],[0 26.66],[32.22 84.73]}; 

oFile = 'C:\Users\crowell\Kahuna\data\sabancaya_5_2018\pulseTrack_analysis\eventAveragedImages_2022-11-02.mat';

gatherImgs = true;
scaleImgs  = false;
fitVz      = true;

%%
if gatherImgs
    for ii = 1:length(avgImgFiles)
        fprintf('%i/%i...\n',ii,length(avgImgFiles))
        I = load(avgImgFiles{ii});
        
        % Load params into stats for image scaling
%         Tavg(ii) = I.Tavg;
        IS.event        = events{ii};
        IS.description  = desc{ii};
        IS.file         = avgImgFiles{ii};
        IS.mask         = I.Tavg.mask;
        IS.maskI        = any(IS.mask,2);
        IS.tSpan        = timeSpans{ii};
        IS.Timg         = I.Tavg.Timg(IS.maskI,:);
        IS.Tfull        = I.Tavg.Tfull(IS.maskI,:);
        IS.Vx           = I.Tavg.Vx;
        IS.VxMasked     = I.Tavg.VxMasked;
        IS.Vz           = I.Tavg.Vz;
        IS.VzMasked     = I.Tavg.VzMasked;
        IS.dTzero       = I.Tavg.dTzero;
        IS.Tz           = max(IS.Timg,[],2); % Go with max for now as a rep core temp
        IS.x            = I.Tavg.x;
        IS.z            = I.Tavg.z(IS.maskI);
        
        IS              = gaussProfileFit(IS,[],[],fitVz);
        
        subT            = IS.Timg;
        subT(~IS.mask(IS.maskI,:)) = NaN;
        IS.prcvals      = I.Tavg.stats.prcvals;
        IS.prctile      = prctile(subT,IS.prcvals,2)';
        IS.mean         = mean(subT,2,'omitnan')';
        IS.var          = nanvar(subT,[],2)';
        IS.max          = nanmax(subT,[],2)';
        IS.min          = nanmin(subT,[],2)';
        IS.saturation   = []; %sum((subT>=tk(ii).satVal),2)';
        IS.nanI         = find(isnan(IS.mean))';
        
        %I.Tavg.Istats = IS;
        Tavg(ii) = IS; %I.Tavg;
    end
    disp('saving...')
    save(oFile,'Tavg')
end

