function yErr = yErrMu(y,w,errMode)
% Function to calculate the average/total uncertainty across track results
%------ y value error function to get mean error-------
% y must be [nx3], with columns equal to [low central high] estimates
% w must be vector of length n. Assumes w = rmse, which is then inverted
%       and normalized
%
% output yErr is [1x3] with average estimates of [low central high] values
%
% (1) Mean and std dev of central estimates
%     yErrMu = @(y) nanmean(y(:,2),1) + [-1 0 1].*nanstd(y(:,2));
%    
% (2) Weighted mean using rmse
%     errWeight = @(x) (1-(x)./max(x(:)));
%     yErrMu = @(y,w) sum(y(:,2).*errWeight(w),1,'omitnan')./sum(errWeight(w),'omitnan') + [-1 0 1].*nanstd(y(:,2),w);
%
% (3) Taking the mean central value and mean magnitude of uncertainty
%       yErrMu = @(y) nanmean(y(:,2)) + [-nanmean(diff(y(:,1:2),[],2),1) 0 nanmean(diff(y(:,2:3),[],2),1)];
%
% (4) Same as (3) but weighted 
%
% () Taking mean central value and mean of absolute low/high bounds
%       yErrMu = @(y) nanmean(y,1);

if nargin<2
    w = [];
end
if nargin<3 && isempty(w)
    errMode = 1;
elseif nargin<3
    errMode = 2;
end

assert(size(y,2)==3,'Check y dims') % Check y dimensions
if ~isempty(w) % Check weight dimensions
    assert(size(w,1)==size(y,1) && size(w,2)==1,'Check w dims')
end

switch errMode
    case 1 % Mean and standard deviation of the track results
        yErr = nanmean(y(:,2),1) + [-1 0 1].*nanstd(y(:,2));
    case 2 % Weighted mean and standard deviation of the track results
%         errWeight = @(x) (1-(x)./max(x(:))); % Normalize the weights
        errWeight = @(x) (1-(x)./(max(x)+min(x))); % Use an inverted normalized rms
        yErr = sum(y(:,2).*errWeight(w),1,'omitnan')./sum(errWeight(w),'omitnan') + [-1 0 1].*std(y(:,2),w,'omitnan');
    case 3
        dy = mean(diff(y,[],2),1);
        yErr = [-dy(1) 0 dy] + nanmean(y(:,2),1);
    case 4
        errWeight = @(x) (1-(x)./(max(x)+min(x))); % Use an inverted normalized rms
        dy = sum(diff(y,[],2) .* errWeight(w),1,'omitnan')./sum(errWeight(w),'omitnan');
        yErr = [-dy(1) 0 dy(2)] + sum(y(:,2).*errWeight(w),1,'omitnan')./sum(errWeight(w),'omitnan');
    case 5 %min max weighted uncertainty for purely random error
        errWeight = @(x) (1-(x)./(max(x)+min(x))); % Use an inverted normalized rms
        [ylo,loi] = nanmin(y,[],1);
        [yhi,hii] = nanmax(y,[],1);
        dy = [(yhi(1)-ylo(1))./2 ...
            (yhi(3)-ylo(3))./2];
        yErr = nanmean(y(:,2),1) + [-dy(1) 0 dy(2)];
    case 6 % Std dev of low and high values
        dy = nanstd(y,0,1);
        yErr = [-dy(1) 0 dy(2)] + nanmean(y(:,2),1);
    case 7 % Weighted mean all
        errWeight = @(x) (1-(x)./(max(x)+min(x)));
        yErr = sum(y.*errWeight(w),1,'omitnan')./sum(errWeight(w),'omitnan');
        
end

end