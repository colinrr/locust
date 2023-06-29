% ---- PROJECT 25B MANUAL DETECTION EDITS -----
loadif(velCube,'V')
loadif(thermCube,'D')

oFile    = fullfile(cubeDir,'v2Tracks/Vtracks_3-6_25B.mat');
% trackFile    = fullfile(cubeDir,'25B_tracks_1_3-6_out_of_date.mat');

figdir = '/Users/crrowell/Kahuna/phd-docs/group-stuff/my-talks/aug-13-2020/';
vidName = 'trackedClusters_25B';
FR = 20;

vidIdx = 1:901;

% 25B Start and end indices to use for each track



runTracker        = true;
saveOutput        = false;
postProcessSmooth = false;
postProcessClip   = false;
plotResults       = true;
clusterVid        = false;

% plotRiseDiagram(D.z,D.t,D.T,'mask',D.mask,'atmo',D.atmo)


%%

    % STA/LTA data input
    thermSource.dataType      = 'var'; % Which field of T0 (source window thermal stats) to use for detection
    thermSource.detectionChannel = 3;
    thermSource.prcvals          = [5 25 50 75 95];

    % Source and detection window parameters
    % thermSource.trackWindowStart = 1;
    % thermSource.trackWindowHeight = 40;
    % thermSource.detectionWindowOffset = 15; % (pixels)
    % thermSource.detectionWindowHeight = 25; % (pixels)
    thermSource.trackWindowStart = 10;
    thermSource.trackWindowHeight = 25;
    thermSource.detectionWindowOffset = 0; % (pixels)
    thermSource.detectionWindowHeight = 25; % (pixels)


    % STA/LTA detection parameters
    detParams.taperLength   = [];
    % detParams.movMinLength  = 15;
    detParams.movMinLength  = 10;
    detParams.l_sta         = [];
    detParams.l_lta         = [];
    detParams.th_on         = 0.44;
    detParams.th_off        = 0.44; 
    detParams.min_dur       = []; 
    detParams.lta_mode      = [];
    detParams.plotflag      = true;
    thermSource.detParams = detParams;

    % Key events and trigger frames
    triggers = {
%        '1' %   1
       '2' %   1
       '3' %   90
%        '4' %  140
%        '5' %  180
%        '6' %  328
%        '7' %  380
    };

if runTracker
    for kk=1:length(triggers)
        clear trackPar
        % These are more efficient to define manually after initial calc
        trackPar.uMax           = 26.82215; 
        trackPar.memoryN        = 9;        
        trackPar.qcPlot         = false;
        
             % MANUAL INPUT SETS FOR INDIVIDUAL 25B pulses
            
        switch triggers{kk}

 
            case '1'
                % Plume front
                trackPar.iTrigger = 1; 
                trackPar.trackWindow = [105 165];
                trackPar.detectWindow = trackPar.trackWindow;
                trackPar.minClust       = 2;
                trackPar.maxClust       = 3; %5
                trackPar.Tpercentile    = 20;
                trackPar.stopTime       = 90;
%                 trackPar.memoryN     = 7;
%                 trackPar.lambda = 0.6;

            case '2'
                % Frame 1, second big cluster - okay performance
                trackPar.iTrigger = 1; 
                trackPar.trackWindow = [41 88]; %[33 72];
                trackPar.detectWindow = trackPar.trackWindow;
                [~,trackPar.prior.mask] = getROI(D.mask(:,:,1),'iLims',trackPar.trackWindow,'jLims',[152 190]);
%                 trackPar.lambda = 0.6;
%                 trackPar.priorWeights   = [0.75 0.25 2.5 0.5];
%                 trackPar.minClust       = 2;
%                 trackPar.maxClust       = 6; %5
                trackPar.stopTime       = 90;

            case '3'
                % 90
                trackPar.iTrigger = 90; 
                trackPar.trackWindow = [1 36]; %[33 72];
                trackPar.detectWindow = trackPar.trackWindow;
                [~,trackPar.prior.mask] = getROI(D.mask(:,:,trackPar.iTrigger),'iLims',[15 36],'jLims',[138 185]);
                trackPar.Tpercentile = 35;
                trackPar.lambda = 0.2;
                trackPar.priorWeights   = [1 0.25 1.75 0.5];
                trackPar.minClust       = 3;
                trackPar.maxClust       = 5; %5
                trackPar.winSzRatio = 1.35;
                trackPar.stopTime = 30;

            case '4'
%                 Excellent ensemble result
                trackPar.iTrigger = 140;
                trackPar.trackWindow = [10 25];
                trackPar.detectWindow = trackPar.trackWindow;
                trackPar.minClust       = 3;
                trackPar.maxClust       = 8; %5
                trackPar.Tpercentile = 50;
                trackPar.memoryN     = 9;
                trackPar.lambda = 0.4;
%                 trackPar.uTol   = 2.75;
                trackPar.clusterWeights = [1.25 0.75 2 1.5 1];
                trackPar.priorWeights = [1.0 0.75 1.75 0.75];
%                 trackPar.winSzRatio = 1.3;
                trackPar.stopTime = 50;

            case '5'
                % 180ish
                trackPar.iTrigger = 275; %180;
                trackPar.trackWindow = [30 75];
                [~,trackPar.prior.mask] = getROI(D.mask(:,:,trackPar.iTrigger),'iLims',[42 71],'jLims',[193 232]);
                trackPar.priorWeights   = [1.5 0.25 1.5 1]; %[1 .25 2 1]
                trackPar.cluterWeights = [0.5 0.75 2 0.5 0.5];
                trackPar.lambda = 0.3;
                trackPar.Tpercentile = 45;
%                 trackPar.uTol = 2.25;
%                 trackPar.minClust = 2;
                trackPar.stopTime = 58;
                
            case '6'
                % 249
                trackPar.iTrigger = 249; 
                trackPar.trackWindow = [1 25];
                [~,trackPar.prior.mask] = getROI(D.mask(:,:,trackPar.iTrigger),'iLims',trackPar.trackWindow,'jLims',[166 193]);
                trackPar.detectWindow = trackPar.trackWindow;
                trackPar.minClust       = 4;
                trackPar.maxClust       = 6;
                trackPar.Tpercentile = 50;
                trackPar.memoryN     = 9;
                trackPar.uTol = 1.75;
%                 trackPar.lambda = 0.5;
                trackPar.stopTime = 70;

            case '7'
                % 380
                trackPar.iTrigger = 380; 
                trackPar.trackWindow = [1 25];
                trackPar.minClust       = 3;
                trackPar.maxClust       = 4;
                trackPar.Tpercentile = 60;
                trackPar.memoryN     = 9;
                trackPar.uTol        = 1.8;
                trackPar.lambda = 0.5;
                trackPar.stopTime = 105;
                
        end
        
    %% DO THE THING
    
        [Vtrack(kk),trackPars(kk)] = trackStructure(V,D,trackPar);
    end

    if saveOutput
        save(oFile,'Vtrack','trackPars')
    end
end
%%

% Smooths structure boundaries between frames, cuts down on high frequency
% noise in statistics retrieval
if postProcessSmooth
    disp('smoothing clusters')
    loadif(trackFile,'Vtrack')
    loadif(trackFile,'trackPars')
    
    for tk=1:length(Vtrack)
        cmask = false([size(D.mask,[1 2]) Vtrack(tk).N]);
        for ti=1:Vtrack(tk).N
            cm = cmask(:,:,ti);
            cm(Vtrack(tk).clustIdx{ti}) = true;
            cmask(:,:,ti) = cm;
        end
        cmask = smooth3(double(cmask),'box',[3 3 D.maskSmoothLength])>0.5;
        for ti=1:Vtrack(tk).N
            Vtrack(tk).clustIdx{ti} = find(cmask(:,:,ti));
        end
    end
end


% This step excludes pixels from occuring in two different structures,
% giving preference to the one occuring later in time
if postProcessClip
    disp('clipping clusters')
    loadif(trackFile,'Vtrack')
    loadif(trackFile,'trackPars')
    
    VtrackPre = Vtrack;
    nT = length(Vtrack);
    
    for tk=1:nT-1
        
        for ti=1:Vtrack(tk).N
            for ci = tk+1:nT
                [tIcheck,tIloc] = ismember(Vtrack(tk).tI(ti),Vtrack(ci).tI);
                if tIcheck
                    Vtrack(tk).clustIdx{ti} = setdiff(Vtrack(tk).clustIdx{ti},Vtrack(ci).clustIdx{tIloc});
                end
            end
            [subI,~] = ind2sub(size(D.T(:,:,1)),Vtrack(tk).clustIdx{ti});
            Vtrack(tk).clustI(ti,:) = [min(subI) mean(subI) max(subI)];
            Vtrack(tk).npx(ti) = numel(Vtrack(tk).clustIdx{ti});
        end
    end
    
end

%% Manual track for comparison
if plotResults


    plotRiseDiagram(D.z,D.t,D.T,'mask',D.mask,'colormap',bone(200),'tracks',Vtrack) %,'atmo',D.atmo)
    caxis([-5 100])

    
%     figure('position',[50 50 1200 1200])
%     plotThermVelocities(D.x,D.z,[],[],'thermal',D.T,'idx',1:10:max(vidIdx),...
%             'mask',D.mask,'tracks',Vtrack,... 'roi',[min(D.z) 1200 -300 700],...
%             'atmo',D.atmo,'Trange',[-5 100],'colormap',gray(200),'time',D.t)
    
end

%% Output cluster video

if clusterVid
    oFile = fullfile(figdir,vidName);
    vidObj = VideoWriter(oFile,'Motion JPEG AVI');
    vidObj.FrameRate = FR;
    open(vidObj);    
    
    fig=figure('position',[50 50 900 1200]);
    for kk=vidIdx
        plotThermVelocities(D.x,D.z,[],[],'thermal',D.T,'idx',kk,...
            'mask',D.mask,'tracks',Vtrack,'roi',[min(D.z) max(D.z) -450 880],...
            'atmo',D.atmo,'Trange',[-5 110],'colormap',gray(200),'time',D.t)

        set(gca,'position',[0.03 0.05 0.95 0.91],'FontSize',12)
        xlabel('X [m from center]')
        zlabel('Z [m above vent]')
        CC=get(gcf,'Children');
        CC(1).Label.FontSize = 12;

        F = getframe(fig);
        writeVideo(vidObj,F);
        clf
    end
    close(vidObj);

end