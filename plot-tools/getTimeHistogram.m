function H = getTimeHistogram(matDir,matHeads,Idx,mask,Tedges)
%   matDir = directory of input .mat frames
%   matHeads = header table asociated with mat files
%
% OPTIONAL input:
%   Idx      = indices to grab. Default: all in matHeads table
%   mask     = optional mask to apply to all frames
%   Tedges   = histogram temperatures. Masked values will be set to the
%               minimum of Tedges.
%               -> empty = min and max Temp taken from first frame
%               -> scalar = minimum Temp. Max Temp taken from first frame
%               -> [1 x 2] = taken as min and max temperaturs, bin size 2
%               -> vector = used directly as bin edges for histcounts.
%
% CROWELL 2020
if nargin<5
    Tedges = [];
end
if nargin<4
    mask = [];
end
if nargin<3
    mask = [];
end
if isempty(mask)
    apply_mask = false;
else
    apply_mask = true;
end

load(matHeads)
if isempty(Idx)
    Idx = str2double(T.Properties.RowNames); % String or double?
end

textprogressbar('Retrieving temperature histograms...')
for kk = 1:length(Idx)
    idx = Idx(kk);
    load(fullfile(matDir,T.File{num2str(idx)}))
    
    if kk==1
        if isempty(Tedges)
            Tedges = round(min(Frame(:))):2:round(max(Frame(:)));
        elseif isscalar(Tedges)
            Tedges = round(Tedges):2:round(max(Frame(:)));
        elseif numel(Tedges) == 2
            Tedges = round(min(Tedges)):2:round(max(Tedges));
        elseif ~isvector(Tedges)
            error('Unrecognized "Tedges" input')
        end
        Tbins    = Tedges(1:end-1) + diff(Tedges);
        nullVal = min(Tedges);
        maxVal  = max(Tedges);
        N = zeros(length(Tbins),length(Idx));
    end
    
    if apply_mask
        Frame(mask) = nullVal;
    end
    N(:,kk) = histcounts(Frame,Tedges);
    textprogressbar(kk/length(Idx)*100)

end
textprogressbar('')

H.Idx = Idx;
H.Tbins = Tbins;
H.Tedges = Tedges;
H.Counts = N;
  
% plotTimeHistogram(H)
end