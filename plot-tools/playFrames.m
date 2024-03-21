function playFrames(matDir,heads,FR,idx,refIdx,roi,cax,maskflag)
% playFrames(matDir,heads,FR,idx,refIdx,roi,cax,maskflag)
% Simple script to load indiviudal .mat file frames and play them video 
% style. Useful for quick scans through a sequence.
%    matDir = directory to .mat frame files
%    heads = path to parameter table containing frame meta-data
% OPTIONAL:
%    FR     = frame rate (per second) [default: 10]
%    idx    = vector of frame indices to play [defaults to all]
%    refIdx = if a reference idx is entered, all frames will be displayed
%             against the reference image using "imshowpair". Leave empty
%             [] to use normal viewing
%    roi    = pixels limits for zooming in [i1 i2 j1 j2]
%    cax    = color axis limits
%    maskflag = true to plot mask outlines if present in header table

if nargin<8
    maskflag = false;
end
if nargin<7
    cax = [];
end
if nargin<6
    roi = [];
end
if nargin<5
    refIdx = [];
end
if nargin<4
    idx = [];
end
if nargin<3
    FR = [];
end

load(heads)
dt = 1./FR;

if ~isempty(refIdx)
    refMode = true;
    load(fullfile(matDir,T.File{num2str(refIdx)}))
    Fref = Frame;
else
    refMode = false;
end

if isempty(idx)
    idx = cellfun(@str2num, T.Properties.RowNames)';
end
load(fullfile(matDir,T.File{idx(1)}))
imsz = size(Frame);
if isempty(roi)
    roi = [1 imsz(1) 1 imsz(2)];
end
if isempty(FR)
    FR = 10;
end


ff=figure('position',[20 20 imsz(2) imsz(1)]);
for ix = idx
    load(fullfile(matDir,T.File{num2str(ix)}))
    ax = subplot(1,1,1);
    if refMode
        imagesc(ax,abs(Fref(roi(1):roi(2),roi(3):roi(4))-Frame(roi(1):roi(2),roi(3):roi(4))))
        colormap(gray(150))
        colorbar
        if ~isempty(cax)
            caxis(cax)
        end
    else
        imagesc(ax,Frame(roi(1):roi(2),roi(3):roi(4)))
        colormap(plasmagrey(150))
        colorbar
        if ~isempty(cax)
            caxis(cax)
        end
    end
    
    if maskflag
        hold on
        if ~exist('mask','var')
            if ismember('Mask',T.Properties.VariableNames)
                mask = T.Mask{num2str(ix)};
            end
        end
        if issparse(mask)
            mask = full(mask);
        end
        if ~islogical(mask)
            mask = logical(mask);
        end
        if any(mask(:))
            mpoly = mask2poly(mask);
            for pp=1:length(mpoly)
                if pp==1;col='w';else; col='c';end
                plot(mpoly(pp).X,mpoly(pp).Y,col,'LineWidth',1.5)
            end
        end
    end
    leg=sprintf('%i\n%.2f s',ix,T.Time(num2str(ix)));
    text(ax,0.2,0.9,leg,'FontSize',12,'Color', [0.9 0.9 0.9],'Units','normalized','HorizontalAlignment','left');

    pause(dt)
    if ix~=idx(end)
        clf
    end
end


end