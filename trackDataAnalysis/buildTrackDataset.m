function [tk,allTrkPar]=buildTrackDataset(tk,thermCubes,velCubes,tracks,errThresh,winLengths,oDir)
%
%
%
%% Test input
% clearvars -except D
% close all

% dataDir   = '~/Kahuna/data/sabancaya_5_2018/';
% tkFile           = fullfile(dataDir,'pulseTrack_analysis/allTracks_2020-10-29_24A_25A4_25B.mat');
% 
% thermCubes = {
%     fullfile(dataDir,'image_exports/24A/thermCubeAnalysis/thermStats_2020-06-29_z644_x578_t1195_revCorr.mat');
%     fullfile(dataDir,'image_exports/25A-4/thermCubeAnalysis/thermStats_2020-07-20_z704_x624_t2439.mat');
%     fullfile(dataDir,'image_exports/25B/thermCubeAnalysis/thermStats_2020-10-26_z710_x571_t1289.mat');
%              };
%          
% tracks = [1 9 20]; %[1 9 20];
% 
% writeOFile = true;
%     oFile = '~/Kahuna/data/sabancaya_5_2018/pulseTrack_analysis/imgTracks_%itracks_%s.mat';

%%
    fprintf('\n========= Reconstruct Track Images and Profiles =========\n')

    if nargin<7
        oDir = '';
        oFile = '';
        writeFlag = false;
    elseif exist(oDir,'dir')
        oFile = strrep(fullfile(oDir,'trackData_%itracks_%s.mat'),'\','\\');
        writeFlag = true;
    else
        warning('Output directory not found.')
        writeFlag = false;
    end
    if nargin<6 || isempty(winLengths)
        winLengths = 1;
    end
    if nargin<5 || isempty(errThresh)
%         errThresh = [D.atmo.zErrThresh D.atmo.zUncThresh 2*(D.atmo.T_halfMax - D.atmo.Tmode)]; 
        errThresh = [];
        zFilter  = false;
    else
        zFilter = true;
    end    
    
    if ischar(tk)
        load(tk)
%         tk = allTracks;
%     elseif isstruct(tk)
%         allTracks = tk;
        assert(and(exist('tk','var') , isstruct(tk)))
    end
    currentCube = '';
    currentVcube = '';

    if isempty(tracks)
        tracks = 1:length(tk);
    end
    
    oFile = sprintf(oFile,length(tracks),datestr(now,'YYYY-mm-dd'));
    %%

    tk = tk(tracks);
%     clear allTracks

    % Cube lists
    tkCube = cell(length(tk),1);
    vkVcube = tkCube;
    for ii=1:length(tk)
        [~,tkCube{ii},~]  = fileparts(tk(ii).thermCube);
        [~,tkVcube{ii},~] = fileparts(tk(ii).velCube);
    end

    fprintf('Processing %i tracks...\n',length(tk))
    for ii=1:length(tk)


        % Load and pre-process new cubes if needed
        [~,currCube,~]  = fileparts(currentCube);
        if or(~exist('D','var'),~strcmp(tkCube{ii},currCube))
            tic
            currentCube = thermCubes{contains(thermCubes,tkCube{ii})};
            currentVcube = velCubes{contains(velCubes,tkVcube{ii})};
            fprintf('Loading and processing data cube(s):\n\t%s\n\t%s\n',currentCube,currentVcube)
            load(currentCube,'D')
            load(currentVcube,'V')
            [~,currCube,~]  = fileparts(currentCube);
           
            
            % Check frame subset
            tracksThisCube = strcmp(tkCube,currCube);
            tkI = sort(unique([tk(tracksThisCube).tI])); % All frames used
            
            % Filter out pixels with large Z Error here
            if zFilter
                disp(' --> Building Height Error Filter...')  

                % Subset for a bit of efficiency
                zErr = zeros(size(D.T));
                zErrMax = zeros([size(D.T,1) 3 size(D.T,3)]);
                [zErr(:,:,tkI),zErrMax(:,:,tkI)] = zErrorEstimation(D.x,D.z+(D.geom.Ztarg-D.geom.Z0),D.mask(:,:,tkI),D.geom,'cylindrical');
                zErrUncertainty = squeeze(zErrMax(:,3,:)-zErrMax(:,1,:)); 

                % Get error threshold flags
                zErrFlags = zErr>errThresh(1);
                zUncFlags = zErrUncertainty>errThresh(2);
                zUncFlags = repmat(permute(zUncFlags,[1 3 2]),[1 size(D.T,2) 1]);

                % Apply error threshold flags
                filtMask  = zErrFlags | zUncFlags; % Remove these pixels
            else
                filtMask = [];
            end
            
            % Find saturated pixels
            disp(' --> Finding saturated pixels...')
            satIdx = cell(size(D.T,3),1);
            for ix = 1:length(satIdx); satIdx{ix} = find(D.T(:,:,ix)>=D.satVal); end
            
            % Remove atmospheric profile from thermal cube
            disp(' --> Removing atmospheric profile...')
            T = zeros(size(D.T));
            mask = false(size(D.mask));
            [T(:,:,tkI),mask(:,:,tkI),dTzero] = removeAtmoProfile(D.mask(:,:,tkI),D.T(:,:,tkI),D.atmo,[],true); %'true');

            if length(winLengths)>1
                winLength = winLengths(contains(thermCubes,tkCube{ii}));
            else
                winLength = winLengths;
            end
            fprintf(' --> Done (%.3f s)\n\n',toc)
        end
        t1 = tic;
        fprintf('Track %i...\n',ii)

        assert(all(tk(ii).imSz==size(D.T,1,2)),'Track and Thermal cube image sizes do not match.')
        

        % Get images with z0, r(z),T(z), W(z) profiles
        tk(ii).dat = getTrackData(T(:,:,tk(ii).tI),mask(:,:,tk(ii).tI),V.Vz(:,:,tk(ii).tI),D.x,D.z,tk(ii),filtMask(:,:,tk(ii).tI),dTzero,satIdx(tk(ii).tI),winLength);  
        fprintf('  --> DONE (%.3f s)\n\n',toc(t1))
        

    end

    %% Write output
    if writeFlag
        fprintf('Writing:\n\t%s\n',oFile)
        save(oFile,'tk','allTrkPar','tracks','errThresh')
    end

end