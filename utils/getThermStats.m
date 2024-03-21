function S = getThermStats(T,M,z,t,subt,subz,prcvals,atmo,satVal,pdfFlag)
% S = getThermStats(T,M,z,t,subt,subz,prcvals,atmo,satVal,pdfFlag)
% Retrieve a profile of statistics from 3D data cube using a 2- or
% 3-D "minicube" sliding window
% INPUT:    T  = data cube
%           M  = data cube mask
%           z  = z vector corresponding to 1st dim of T and M
%           t  = time vector corresponding to 3rd dim of T and M
%           subt = [N samples x P windows] TIME indices for each window
%           subz = array: [M samples x P windows] HEIGHT indices for each window
%               OR cell: length Px1, with each entry containing N pixel 
%                 indices for the corresponding frame
% OPTIONAL IN:
%           prcvals = vector of percentile values to calculate. 
%                     Default: [5 25 50 75 95]
%           atmo = Struct output from fitAtmoProfile. Remove atmospheric
%                   profile and apply Temperature cutoff to mask, with 
%                   default nHalfMax = 2, unless nHalfMax is specified as a
%                   field in 'atmo'.
%           satVal = Saturation temperature value. Data with this value
%                   will be flagged. 
%
% OUT fields: mean, median, variance/std, percentiles, max, min for each window 
%       histcounts?
%
%  C Rowell Sep 2019
% NOTES: Implement later: Input:  Tbins = [optional] bin edges to pull T histcounts

nHalfMax = 2;

if nargin<10
    pdfFlag = false;
end
if nargin<9
    satVal = [];
end
if nargin<8
    atmo = [];
end
if nargin<7
    prcvals = [];
end

if isempty(prcvals)
    prcvals = [5 25 50 75 95];
end
if isempty(atmo)
    removeAtmo = false;
else
    removeAtmo = true;
end

assert(numel(t)==size(T,3),'t vector must match 3rd dimension of T and M')
assert(all(size(T)==size(M)),'T and M must have the same size')

if and(size(t,1)>1,isvector(t))
    t = t';
end
if and(size(z,1)>1,isvector(z))
    z = z';
end

NumWins = size(subt,2);

if iscell(subz)
    assert(size(subt,2)==length(subz),'subt and subz must have the same length/dimension corresponding to number of frames.')
    idxMode = true;
    NumVals = [];
    Tcut = cell(NumWins,1);

elseif isnumeric(subz)
    assert(size(subt,2)==size(subz,2),'subt and subz must have the same length in dimension two, corresponding to number of frames.')
    idxMode = false;
    NumVals = size(subt,1).*size(subz,1).*size(T,2); % Number of values per mini cube
    Tcut = nan(NumVals,NumWins);
end

if ~isempty(satVal)
    satM = T(:,:,subt)>=satVal; % Saturation mask
    satFlag = true;
else
    satFlag = false;
end

if removeAtmo % Profile removal
    if isfield(atmo,'nHalfMax')
        nHalfMax = atmo.nHalfMax;
    end
    [T,M] = removeAtmoProfile(M,T,atmo,nHalfMax,true);
%     T = T - atmo.Tinterp - atmo.T_halfMax; % Profile removal
%     M = and(M, T>= -nHalfMax*(atmo.Tmax - atmo.T_halfMax)); % Mask threshold
end

% Be nice to vectorize this...

switch idxMode
    % USING INDICES
    case true
        S.tI        = subt;
        S.zI        = subz;
        S.t         = mean(t(subt),1);
        S.z         = zeros(1,NumWins);
        S.N         = zeros(1,NumWins);
        S.prcvals   = prcvals;
        S.prctile   = zeros(length(prcvals),NumWins);
        S.mean      = zeros(1,NumWins);
        S.var       = zeros(1,NumWins);
        S.max       = zeros(1,NumWins);
        S.min       = zeros(1,NumWins);
        if satFlag
            S.saturation = zeros(1,NumWins);
        else
            S.saturation = [];
        end
        
        for kk=1:NumWins    
            Tkk = T(:,:,subt(kk));
            mkk = M(:,:,subt(kk));
            % Remove wherever there is no mask (no detected object)
            idx = intersect(find(mkk),subz{kk});
            S.N(kk) = numel(idx);
            Tcut{kk} = Tkk(idx);            
            if satFlag
                satMkk = satM(:,:,kk);
                S.saturation(kk) = sum(satMkk(idx));
            end
            
%             if removeAtmo
% %                 mkk(Tkk< (-nHalfMax*(atmo.Tmax - atmo.T_halfMax)) ) = false;
%                 idx = intersect(find(mkk),subz{kk});
%                 Tcut{kk} = Tkk(idx);
%             else
% 
% 
%             end
%             Tcut{kk}(Tcut{kk}<(-nHalfMax*(atmo.Tmax - atmo.T_halfMax))) = NaN;
        
            [subI,~]     = ind2sub(size(Tkk),subz{kk});
            S.z(kk)         = mean(z(subI));
            S.prctile(:,kk) = prctile(Tcut{kk},S.prcvals);
            S.mean(kk)      = nanmean(Tcut{kk});
            S.var(kk)       = nanvar(Tcut{kk});
            S.max(kk)       = nanmax(Tcut{kk});
            S.min(kk)       = nanmin(Tcut{kk});
        
        end
        
    % USING Z SUBSCRIPTS    
    case false
        satVec = zeros(1,NumWins);
        for kk=1:NumWins    
            subzk = subz(~isnan(subz(:,kk)),kk);
            Tminicube = T(subzk,:,subt(:,kk));
            % Place NaN's wherever there is no mask (no detected object)
            Tminicube(~M(subzk,:,subt(:,kk))) = NaN;
            Tcut(1:numel(Tminicube),kk) = Tminicube(:);
            if satFlag
                satMkk = satM(subzk,:,kk);
                satVec(kk) = sum(satMkk(M(subzk,:,subt(:,kk))));
            end
            
            % Crude way to get around window nans
            if ~any(isnan(subz(:,kk)))
                S.z(kk) = mean(z(subz(:,kk)));
            else
                nn = ~isnan(subz(:,kk));
                ix = 1:size(subz,1);
                S.z(kk) = interp1(ix(nn),z(subz(nn,kk)),mean(ix),'linear','extrap');
            end
        end
        
        S.tI        = subt;
        S.zI        = subz;
        S.t         = mean(t(subt),1);
%         S.z         = mean(z(subz),1);
        S.N         = sum(~isnan(Tcut),1);
        S.prcvals   = prcvals;
        S.prctile   = prctile(Tcut,S.prcvals,1);
        S.mean      = nanmean(Tcut,1);
        S.var       = nanvar(Tcut,[],1);
        S.max       = nanmax(Tcut,[],1);
        S.min       = nanmin(Tcut,[],1);
        if satFlag
            S.saturation = satVec;
        else
            S.saturation = [];
        end
end

% for kk=1:NumWins
%     if idxMode
% %         [subI,subJ] = ind2sub(size(T(:,:,1),[1 2]),subz{kk});
%         Tkk = T(:,:,kk);
%         Tcut{kk} = Tkk(subz{kk});
%     else
%         Tminicube = T(subz(:,kk),:,subt(:,kk)).*M(subz(:,kk),:,subt(:,kk));
%         Tcut(:,kk) = Tminicube(:);
%     end
% end

% Populate output struct

% if idxMode

%     for  jj=1:NumWins
%         Tcut{jj}(Tcut{jj}<(-nHalfMax*(atmo.Tmax - atmo.T_halfMax))) = NaN;
        
%     end
% else
    % Place NaN's wherever there is no mask (no detected object)
%     Tcut(Tcut==0) = NaN;
% 
%     S.z         = mean(z(subz),1);
%     S.prcvals   = prcvals;
%     S.prctile   = prctile(Tcut,S.prcvals,1);
%     S.mean      = nanmean(Tcut,1);
%     S.var       = nanvar(Tcut,[],1);
%     S.max       = nanmax(Tcut,[],1);
%     S.min       = nanmin(Tcut,[],1);
    % S.subz      = subz;
    % S.subt      = subt;
% end

% Follow up NaN check
N1 = isnan(S.prctile);
N2 = isnan(S.mean);
N3 = isnan(S.var);
N4 = isnan(S.max);
N5 = isnan(S.min);

tf = isequaln(N1(1,:),N2,N3,N4,N5);
if ~tf
    warning('NaN values should only present where there is no mask, but NaN locations are not equal between fields.')
end
S.nanI = find(N2);

% S.prctile()   = 0;
% S.mean(isnan(S.mean))      = 0;
% S.var(isnan(S.var))       = 0;
% S.max(isnan(S.max))       = 0;
% S.min(isnan(S.min))       = 0;


end