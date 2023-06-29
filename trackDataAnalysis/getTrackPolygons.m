function tk = getTrackPolygons(tk,eventMasks,oFile,plotFlag) %,maskFile)
% Attaches polygons to track file - frame by frame for plume and cluster
% masks, plus a single boundary polygon for the entire track over time
% tk         = track struct or path to struct file
% eventMasks = struct of plume masks or path to file
% oFile      = output path for writing file. Only works if tk is also a path

if nargin<4
    plotFlag = false;
end
if nargin<3
    oFile = [];
end
    

if isempty(oFile)
    writeFlag = false;
elseif ~isempty(oFile)  && ischar(tk) && exist(tk,'file')
    load(tk)
    writeFlag = true;
end

if ischar(eventMasks)
    load(eventMasks)
end


maskThresh = 1;

textprogressbar('   Getting track polygons:  ')
textprogressbar(0)

for ti = 1:length(tk)
    
    nf = length(tk(ti).tI); % Num frames in track
    ei = find(ismember({eventMasks.name},tk(ti).event)); % Which event mask?
    cmask0 = false(eventMasks(ei).size([1 2]));
    cmaskAll = false([eventMasks(ei).size([1 2]) nf]);
    
    for ii = nf:-1:1
        [~,~,pPoly(ii)] = getROI(eventMasks(ei).M(:,:,tk(ti).tI(ii)),'maxRegions',1);
        cmask = cmask0;
        cmask(tk(ti).clustIdx{ii}) = true;
        cmaskAll(:,:,ii) = cmask;
        [~,~,tPoly(ii)] = getROI(cmask,'maxRegions',1);
    end
    
    maskStack = sum(double(cmaskAll),3);
    [trackROI,trackMask,trackBounds] = getROI(maskStack>=maskThresh,'maxRegions',1);
    
    tk(ti).x = eventMasks(ei).x;
    tk(ti).z = eventMasks(ei).z;
    tk(ti).x0 = eventMasks(ei).x0;
    tk(ti).z0 = eventMasks(ei).z0;
    tk(ti).trackPoly = tPoly;
    tk(ti).plumePoly = pPoly;
    tk(ti).trackBounds = trackBounds;
    tk(ti).trackBounds.mask = trackMask; 
    tk(ti).trackBounds.roi  = trackROI;
    
    textprogressbar(ti/length(tk)*100)
    clear trackBounds tPoly pPoly
end
textprogressbar(' --> Done')
if plotFlag
    nr = floor(sqrt(length(tk)));
    nc = ceil(sqrt(length(tk)))+1;
    if (nr*nc)<length(tk); nc = nc+1;end
    figure('position',[50 50 1000 800])
    for ti = 1:length(tk)
        subplot(nr,nc,ti)
        plotTrackOutlines(tk(ti))
        title(sprintf('%i: %s-%i',ti,tk(ti).event,tk(ti).eventTrack))
    end
    
end

if writeFlag
    fprintf('Writing:\n\t%s\n',oFile)
    if and( exist('allTrkPar','var'), exist('errThresh','var') )
        save(oFile,'tk','errThresh','allTrkPar')
    else
        save(oFile,'tk')
    end
end

end