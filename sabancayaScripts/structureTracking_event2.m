% ---- PROJECT 24A MANUAL DETECTION/TRACKING -----
loadif(velCube,'V')
loadif(thermCube,'D')

figdir = '/Users/crrowell/Kahuna/phd-docs/group-stuff/my-talks/aug-13-2020/';
ofile = fullfile(cubeDir,'24A_tracks_1-6'); % Track output

vidName = 'trackedClusters_24A_smoothed'; % Video output
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
% thermSource.prcvals          = [5 20 50 80 95];
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
   '1' %   9
   '2' %   142
   '3' %   205
   '4' %  252
   '5' %  286
%    '6a' %  328
   '6b' %  385
   '6c' %  401
   '7' %   492
% %    '8' %   556
};

% For each case, if there are two sets of parameters, the first is for the
% original tracks with trackVelocities V1.10, the second are updated for
% V1.11
if runTracker
    for kk=1:length(triggers)
        clear trackPar
        trackPar.uMax           = 13.32; % More efficient
        trackPar.memoryN        = 7;        % More efficient

        switch triggers{kk}
            case '1'
    %             Plume front attempt 1 - good
                trackPar.iTrigger = 55;
                trackPar.trackWindow = [1 40];
                trackPar.minClust       = 2;
                trackPar.maxClust       = 3;
                trackPar.lambda          = 0.1;
                trackPar.winSzRatio     = 1.7;
                trackPar.memoryN     = 7;
                trackPar.Tpercentile = 30;
                trackPar.clusterWeights = [1 1 2   0   .25];
                trackPar.priorWeights   = [1 1 1 1];
                trackPar.uTol           = 3;
                trackPar.stopTime       = 175;

            case '2'
                % Pulse 2 - good
                trackPar.iTrigger = 142;
                trackPar.trackWindow = [1 30];
                [~,trackPar.prior.mask] = getROI(D.mask(:,:,trackPar.iTrigger),'iLims',trackPar.trackWindow,'jLims',[95 120]);
                % trackPar.memoryN        = 7;
                trackPar.minClust       = 2;
                trackPar.maxClust       = 4;
                % trackPar.winSzRatio     = 3;
                trackPar.Tpercentile    = 30;
                trackPar.lambda         = 0.4;
                trackPar.uTol           = 4;
%                 trackPar.priorWeights   = [2 2 1 1];
                trackPar.stopTime       = 51;

            case '3'
%                 % Pulse 3 - good
%                 trackPar.iTrigger = 205;
%                 trackPar.trackWindow = [1 20];
%                 % trackPar.winSzRatio     = 1.5;
%                 trackPar.minClust       = 2;
%                 trackPar.maxClust       = 5;
%                 trackPar.stopTime       = 80;

                % Pulse 3 - good
                trackPar.iTrigger = 205;
                trackPar.trackWindow = [1 20];
                [~,trackPar.prior.mask] = getROI(D.mask(:,:,trackPar.iTrigger),'iLims',trackPar.trackWindow,'jLims',[74 94]);
                trackPar.Tpercentile    = 35;
                % trackPar.winSzRatio     = 1.5;
                trackPar.minClust       = 2;
                trackPar.maxClust       = 4;
                trackPar.clusterWeights = [0.5 0.5 1.5 0.5 1.5];
%                 trackPar.priorWeights   = [1 0.25 2 0.5];
                trackPar.uTol           = 3;
                trackPar.winSzRatio     = 1.6;
                trackPar.stopTime       = 80;

            case '4'
    %             Pulse 4a - dope lil' un
                trackPar.iTrigger = 252;
                trackPar.trackWindow = [1 20];
                trackPar.minClust       = 3;
                trackPar.maxClust       = 5;
                trackPar.lambda         = [0.6];
                trackPar.stopTime       = 80;

            case '5'
%    %             Pulse 4b - good
%                 trackPar.iTrigger = 285;
%                 trackPar.trackWindow = [1 20];
%                 trackPar.minClust       = 3;
%                 trackPar.maxClust       = 5;
%                 trackPar.clusterWeights = [1 1 2   0.5   2];
%                 trackPar.lambda         = [0.6];
%                 trackPar.priorWeights   = [2 2 1 1];
%                 trackPar.stopTime       = 100;
                
                %  Pulse 5 - proper good
                trackPar.iTrigger = 290;
                trackPar.trackWindow = [1 20];
                trackPar.minClust       = 2;
                trackPar.maxClust       = 5;
                trackPar.Tpercentile    = 35;
                trackPar.clusterWeights = [1 1 2   0.5   1];
                trackPar.lambda         = [0.6];
%                 trackPar.priorWeights   = [2 2 1 1];
                trackPar.stopTime       = 100;

            case '6a'
                % Pulse 5a
%                 trackPar.iTrigger = 328;
%                 trackPar.trackWindow = [1 10];
%                 trackPar.minClust       = 4;
%                 trackPar.uTol           = 2;
%                 trackPar.stopTime       = 110;
                
                % Pulse 5a
                trackPar.iTrigger = 335;
                trackPar.trackWindow = [1 10];
                trackPar.Tpercentile    = 35;                
                trackPar.minClust       = 2;
                trackPar.maxClust       = 5;
                trackPar.clusterWeights = [0.75 0.75 1   0.5   1];
%                 trackPar.priorWeights   = [0.5 0.25 2 0.5];
                trackPar.lambda         = 0.4;
%                 trackPar.uTol           = 2;
                trackPar.stopTime       = 110;

            case '6b'
%                 % Pulse 5b - loses it  t~108/110?
%                 trackPar.iTrigger = 358;
%                 trackPar.trackWindow = [1 10];
%                 trackPar.minClust       = 4;
%                 trackPar.maxClust       = 7;
%                 trackPar.lambda         = [0.6];
%                 trackPar.uTol           = 2;
%                 trackPar.priorWeights   = [2 0.5 2 1];
%                 trackPar.stopTime       = 115;

                % Pulse 5b - loses it  t~108/110?
                trackPar.iTrigger = 358;
                trackPar.trackWindow = [1 10];
                trackPar.Tpercentile    = 35;                
                trackPar.minClust       = 2;
                trackPar.maxClust       = 5;
                trackPar.lambda         = [0.6];
                trackPar.uTol           = 2;
                trackPar.priorWeights   = [1 0.5 2 1];
                trackPar.stopTime       = 115;

            case '6c'
                % Pulse 5c
%                 trackPar.iTrigger = 385;
%                 trackPar.trackWindow = [1 10];
%     %             trackPar.minClust       = 4;
%     %             trackPar.maxClust       = 7;
%                 trackPar.lambda         = [0.6];
%                 trackPar.uTol           = 2;
%                 trackPar.priorWeights   = [2 0.5 2 1];
%                 trackPar.stopTime       = 108;

                % Pulse 5c
                trackPar.iTrigger = 385;
                trackPar.trackWindow = [1 10];
                trackPar.Tpercentile    = 35;                
                trackPar.minClust       = 2;
                trackPar.maxClust       = 3;
                trackPar.lambda         = [0.3];
                trackPar.clusterWeights = [.75 1 1.   0.25   0.75];
                trackPar.uTol           = 2.2;
                trackPar.priorWeights   = [0.75 0.5 2 0.5];
                trackPar.stopTime       = 120;
                %5d 399 - short-lived bit...

            case '7' % good
                trackPar.iTrigger = 492;
                trackPar.trackWindow = [1 20];
                trackPar.minClust       = 4;
                trackPar.uTol           = 2;
                trackPar.stopTime       = 142;


            case '8' % eh iffy target
    %             trackPar.iTrigger = 399;
    %             trackPar.trackWindow = [1 20];
    %             trackPar.minClust       = 4;
    %             trackPar.stopTime       = 120;

        end
    %% DO THE THING


        [Vtrack(kk),trackPars(kk)] = pulseTrack(V,D,trackPar);
    end

    save(rawTrackFile,'Vtrack','trackPars')
end

%%
if postProcessSmooth
    disp('smoothing clusters')
    loadif(ofile,'Vtrack')
    loadif(ofile,'trackPars')
    
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

if postProcessClip
    disp('clipping clusters')
    loadif(ofile,'Vtrack')
    loadif(ofile,'trackPars')
    
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
    loadif(ofile,'Vtrack')
    loadif(ofile,'trackPars')
    
    plotRiseDiagram(D.z,D.t,D.T,'mask',D.mask,'atmo',D.atmo,'tracks',Vtrack)
    colormap(gray(200))
    caxis([0 70])
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




