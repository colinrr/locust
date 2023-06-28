% ---- PROJECT 25A4 MANUAL DETECTION/TRACKING -----
loadif(velCube,'V')
loadif(thermCube,'D')

% figdir = '/Users/crrowell/Kahuna/phd-docs/group-stuff/my-talks/aug-13-2020/';


vidName = 'trackedClusters_25A4_smoothed'; % Video output
FR = 10;

vidIdx = 1:640;


runTracker        = true;
postProcessSmooth = false;
postProcessClip   = false;
plotResults       = true;
clusterVid        = false;
%%
% STA/LTA data input
% thermSource.dataType      = 'var'; % Which field of T0 (source window thermal stats) to use for detection
thermSource.dataType      = 'prctile'; 
thermSource.detectionChannel = 5;
thermSource.prcvals          = [5 25 50 75 95];

% Source and detection window parameters
thermSource.trackWindowStart = 1;
thermSource.trackWindowHeight = 20;
thermSource.detectionWindowOffset = 0; % (pixels)
thermSource.detectionWindowHeight = 20; % (pixels)
thermSource.removeAtmo = true;

% STA/LTA detection parameters
detParams.taperLength   = [];
% detParams.movMinLength  = 15;
detParams.movMinLength  = 10;
detParams.l_sta         = [];
detParams.l_lta         = [];
detParams.th_on         = 0.7;
detParams.th_off        = []; 
detParams.min_dur       = []; 
detParams.lta_mode      = [];
detParams.plotflag      = true;
thermSource.detParams = detParams;



% Key events and trigger frames
triggers = {
%    '1' %   330
%    '2' %   
   '3' %   
%    '4' %  
%    '5' %  
%    '6' %  
%    '7' %  
%    '8' %  
%    '9' %  
%    '10' %   
%    '11'
};

if runTracker
    for kk=1:length(triggers)
        clear trackPar
        trackPar.uMax           = 8.59; % More efficient
        trackPar.memoryN        = 5;        % More efficient

        switch triggers{kk}
            case '1'
                % Pretty good
                trackPar.iTrigger       = 330;
                trackPar.trackWindow    = [1 20];
                trackPar.minClust       = 2;
                trackPar.maxClust       = 4;
                trackPar.lambda         = 0.2;
%                 trackPar.winSzRatio     = 1.7;
%                 trackPar.memoryN     = 7;
                trackPar.Tpercentile    = 50;
                trackPar.clusterWeights = [1 1 1 0.25 .5];
                trackPar.priorWeights   = [1 .25 1.5 1];
                trackPar.uTol           = 4;
                trackPar.nHalfMax       = 3;
                trackPar.stopTime       = 110;

            case '2'
                % Too small earlier, but ok until t=86.23
                trackPar.iTrigger       = 445;
                trackPar.trackWindow    = [1 20];
                trackPar.minClust       = 2;
                trackPar.maxClust       = 4;
                trackPar.lambda         = 0.4;
                trackPar.Tpercentile    = 40;
                trackPar.clusterWeights = [0.75 1 1.5 0.25 1];
                trackPar.priorWeights   = [1.25 .5 1.5 1];
                trackPar.winSzRatio     = 1.4;
                trackPar.uTol           = 3.3;
                trackPar.nHalfMax       = 3;
                trackPar.stopTime       = 125;

            case '3'
                % Good - nice long track
                trackPar.iTrigger       = 570;
                trackPar.trackWindow    = [1 20];
                trackPar.minClust       = 2;
                trackPar.maxClust       = 5;
                trackPar.lambda         = 0.4;
                trackPar.Tpercentile    = 40;
                trackPar.clusterWeights = [1 1 1.75 0.25 1];
                trackPar.priorWeights   = [1.75 .5 1.5 1.25];
                trackPar.uTol           = 4.5;
                trackPar.nHalfMax       = 3;
                trackPar.stopTime       = 180;
                

            case '4'
                % Not bad...incorporates following pulse but can be
                % excluded
                trackPar.iTrigger       = 680;
                trackPar.trackWindow    = [1 20];
                trackPar.minClust       = 2;
                trackPar.maxClust       = 3;
                trackPar.lambda         = 0.2;
                trackPar.Tpercentile    = 35;
                trackPar.clusterWeights = [0.6 0.6 2 0.5 1.5];
%                 trackPar.priorWeights   = [1.5 .5 1.5 1.25];
                trackPar.uTol           = 4.;
                trackPar.nHalfMax       = 3;
                trackPar.stopTime       = 196;

            case '5'
                % Good early, too small late
                trackPar.iTrigger       = 990;
                trackPar.trackWindow    = [1 20];
                trackPar.minClust       = 2;
                trackPar.maxClust       = 3;
                trackPar.lambda         = 0.2;
                trackPar.Tpercentile    = 40;
                trackPar.uTol           = 3.2;
                trackPar.clusterWeights = [0.6 0.6 2 1 1.];
                trackPar.priorWeights   = [1. .1 1. 1.25];
                trackPar.winSzRatio     = 1.3;
                trackPar.nHalfMax       = 3;
                trackPar.stopTime       = 180;
           
            case '6'
                % Ok in middle, not great early and late
                trackPar.iTrigger       = 1210;
                trackPar.trackWindow    = [1 20];
                trackPar.minClust       = 2;
                trackPar.maxClust       = 3;
                trackPar.lambda         = 0.2;
                trackPar.Tpercentile    = 30;
                trackPar.uTol           = 3.5;
                trackPar.clusterWeights = [0.6 0.6 1.5   0.25   1.5];
%                 trackPar.priorWeights   = [1. .5 1. 1.25];
                trackPar.nHalfMax       = 3;
                trackPar.stopTime       = 250;                

            case '7'
                % Good, but needs trimming from next pulse
                trackPar.iTrigger       = 1346;
                trackPar.trackWindow    = [1 20];
                trackPar.minClust       = 2;
                trackPar.maxClust       = 3;
                trackPar.lambda         = 0.1;
                trackPar.Tpercentile    = 45;
                trackPar.uTol           = 3.5;
                trackPar.clusterWeights = [0.75 0.75 2  1  1.5];
                trackPar.priorWeights   = [1 0.25 2 0.5];
                trackPar.nHalfMax       = 3;
                trackPar.stopTime       = 250;                

            case '8'
                % Hard to pin down, so used combined results of 2 runs. 
                % See 25A5_track8_combined.mat and sandbox_scripts/Track8_25A4.m
                trackPar.iTrigger       = 1495;
                trackPar.trackWindow    = [1 20];
                [~,trackPar.prior.mask] = getROI(D.mask(:,:,trackPar.iTrigger),'iLims',trackPar.trackWindow,'jLims',[40 65]);
                trackPar.minClust       = 2;
                trackPar.maxClust       = 3;
                trackPar.lambda         = 0.2;
                trackPar.Tpercentile    = 42;
                trackPar.uTol           = 3.5;
                trackPar.clusterWeights = [0.75 1. 2  0.75  1];
                trackPar.nHalfMax       = 3;
                trackPar.winSzRatio     = 1.3;
                trackPar.stopTime       = 250;                

            case '9'
                % Good early-ish to mid, then loses diameter. T stays great late
                trackPar.iTrigger       = 1695;
                trackPar.trackWindow    = [1 25];
                [~,trackPar.prior.mask] = getROI(D.mask(:,:,trackPar.iTrigger),'iLims',trackPar.trackWindow,'jLims',[40 77]);
                trackPar.minClust       = 2;
                trackPar.maxClust       = 3;
                trackPar.lambda         = 0.15;
                trackPar.Tpercentile    = 35;
                trackPar.uTol           = 3.5;
                trackPar.clusterWeights = [0.5 0.75 1.5 0.25 1.5];
                trackPar.priorWeights   = [1 0.25 1.5 0.5];
                trackPar.nHalfMax       = 3;
%                 trackPar.stopTime       = 250;  

            case '10'
                % 
                trackPar.iTrigger       = 1828;
                trackPar.trackWindow    = [1 20];
                trackPar.minClust       = 2;
                trackPar.maxClust       = 3;
                trackPar.lambda         = 0.15;
                trackPar.Tpercentile    = 35;
                trackPar.uTol           = 3.5;
                trackPar.clusterWeights = [0.5 0.75 1.5 0.25 1.5];
                trackPar.priorWeights   = [1 0.25 1.5 0.5];
                trackPar.nHalfMax       = 3;

            case '11'
                % 
                trackPar.iTrigger       = 2105;
                trackPar.trackWindow    = [1 20];
%                 [~,trackPar.prior.mask] = getROI(D.mask(:,:,trackPar.iTrigger),'iLims',trackPar.trackWindow,'jLims',[40 77]);
                trackPar.minClust       = 2;
                trackPar.maxClust       = 3;
                trackPar.lambda         = 0.1;
                trackPar.Tpercentile    = 40;
                trackPar.uTol           = 4.0;
                trackPar.clusterWeights = [0.6 0.6 1. 0.25 1.5];
                trackPar.priorWeights   = [1 0.25 1.5 0.5];
                trackPar.winSzRatio     = 1.4;
                trackPar.nHalfMax       = 3;

        end
    %% DO THE THING


        [Vtrack(kk),trackPars(kk)] = pulseTrack(V,D,trackPar);
    end

    save(rawTrackFile,'Vtrack','trackPars')
end

%%
% Smooths structure boundaries between frames, cuts down on high frequency
% noise in statistics retrieval
if postProcessSmooth
    disp('smoothing clusters')
    loadif(rawTrackFile,'Vtrack')
    loadif(rawTrackFile,'trackPars')
    
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
    loadif(rawTrackFile,'Vtrack')
    loadif(rawTrackFile,'trackPars')
    
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
%%

if plotResults
    loadif(rawTrackFile,'Vtrack')
    loadif(rawTrackFile,'trackPars')
    
    plotRiseDiagram(D.z,D.t,D.T,'mask',D.mask,'atmo',D.atmo,'tracks',Vtrack,'idx',1:2:size(D.T,3),'cax',[0 40])
    colormap(gray(200))
    ylim([min(D.z) 1000])
%     caxis([0 70])
end


%%

if clusterVid
    oFile = fullfile(figdir,vidName);
    vidObj = VideoWriter(oFile,'Motion JPEG AVI');
    vidObj.FrameRate = FR;
    open(vidObj);    
    
    fig=figure('position',[50 50 1100 1200]);
    for kk=vidIdx
        plotThermVelocities(D.x,D.z,[],[],'thermal',D.T,'idx',kk,...
            'mask',D.mask,'tracks',Vtrack,'roi',[min(D.z) 1200 -300 700],...
            'atmo',D.atmo,'Trange',[-5 75],'colormap',gray(200),'time',D.t)

        set(gca,'position',[0.03 0.05 0.95 0.90],'FontSize',12)
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




