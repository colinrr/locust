function varargout = ksdensityND(A,pts,dim)
% [f,xi,fAll,stats] = ksdensityND(A,pts,dim)
% Get probability density function for a 1,2, or 3D matrix along one or two
% dimensions. As per the ksdensity default, NaN's are ignored. This can be
% very slow when operating along only a single dimension for a large 3D
% matrix.
%
% IN:   A   = 2D or 3D array
%       pts = specifies points to evaluate f. SAME AS for ksdensity.
%           Here, xi and pts contain identical values.
%       dim = dimension(s) along which to operate in A. Default = 1.      
%
% OUT:  f    = 2 or 3D matrix, where size(f,dim) = length(pts).
%       xi   = pts, as in ksdensity
%       fAll = ksdensity function for all points in A (ie ksdensity(A(:)))
%           E.g. for a 3D array:
%            1) dim = 1 produces a 3D array with sizes [length(xi) size(A,2) size(A,3)]
%                ie. each column of f along dimension 1 is now the probability
%                distrubution function of the corresponding column in A.
%           
%            2) dim = [2 3] produces a 2D array with sizes [size(A,1) length(xi)]
%               ie. each row of f is the probaility distribution of the
%               corresponding 2D slice (dimensions 2 and 3) of A.
%       stats  = struct containing a range of common statistics for each
%            PDF in f
%           Each field is a matrix with values corresponding to the PDF's
%           in f. Output values are: N (number of non-NaN data points),
%           Mean, Median, Mode, Standard Deviation, Min, Max.
%
%
% C Rowell Jun 2020

% Procudes a 2D matrix where the columns (or rows) are the probability density
% functions along dimension 1 (or 2, respectively).

if nargin<3
    dim = 1;
end
if nargin<2
    pts = [];
end

if nargout>=3
    if isempty(pts)
        [varargout{3},pts] = ksdensity(A(:));
    else
        varargout{3} = ksdensity(A(:),pts);
    end
elseif isempty(pts)
    pts = linspace(min(A(:)),max(A(:)),100);
end

nxi = length(pts);
nDimsA = sum(size(A)>1);
allDimsA = (1:nDimsA)';
ndims  = length(dim);
Adims = size(A);
assert(ndims<nDimsA,'For dim along all dimensions of A, just use ksdensity(A(:)).')

% Permute,reshape A to prep
nonDims = allDimsA(~any(allDimsA == dim,2))';
A = reshape(permute(A,[dim nonDims]),[prod(size(A,dim)) prod(size(A,nonDims))]);

% Loop ksdensity
f = zeros(nxi,size(A,2));
% N = zeros(size(A,2),1);
for ii=1:size(A,2)
    if sum(~isnan(A(:,ii)))>0
        f(:,ii) = ksdensity(A(:,ii),pts);
    else
        f(:,ii) = 0;
    end
end

% f = reshape(f,[nxi Adims(nonDims)]);
% Reshape back to fit original dimensions
[~,dI]=sort([dim nonDims]);
[~,maxI] = max(f,[],1);
fmode    = pts(maxI);
varargout{1} = squeeze(permute(reshape(f,[nxi Adims(nonDims)]),dI));
% N = reshape(N,Adims(nonDims));



if nargout>=2
    varargout{2} = pts;
end
if nargout==4
    if length(nonDims)==1
        statDims = [Adims(nonDims) 1];
    elseif length(nonDims)==1
        statDims = Adims(nonDims);
    end
    stats.N      = reshape(sum(A,1,'omitnan'),statDims);
    stats.Mean   = reshape(mean(A,1,'omitnan'),statDims);
    stats.Median = reshape(median(A,1,'omitnan'),statDims);
    stats.Mode   = reshape(fmode,statDims);
    stats.StdDev = reshape(std(A,[],1,'omitnan'),statDims);
    stats.prc25_75 = reshape(prctile(A,[25 75],1),[statDims 2]);
    stats.Min    = reshape(min(A,[],1,'omitnan'),statDims);
    stats.Max    = reshape(max(A,[],1,'omitnan'),statDims);
    
    varargout{4} = stats;
end


end