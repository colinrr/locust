function rStats = getRparams(z,r,idx)
% rStats = getRparams(z,r,idx)
% Fitting function to estimate feature spreading rate and virtual origin
% using a linear fit on radius.
% INPUT: z and radius vectors. idx = optional vector subset of indices

    if nargin<3
        idx = 1:length(z);
    end
    if numel(idx)==2
        idx = idx(1):idx(2);
    end
    
    rStats.zLim = [z(idx(1)) z(idx(end))];
    rStats.z = z(idx);
    rStats.r = r(idx);
    rStats.idx = idx;
    
    % Linear model for R = f(z)
    rStats.Rmdl = fitlm(rStats.z,rStats.r,'linear','RobustOpts','on');
    rStats.r0   = rStats.Rmdl.Coefficients.Estimate(1);
    coefs       = coefCI(rStats.Rmdl);
    rStats.r0ci = coefs(1,:);
    rStats.drdz = rStats.Rmdl.Coefficients.Estimate(2);
    rStats.drdzCI = coefs(2,:);
    rStats.grad_r_z = gradient(rStats.r,rStats.z);
    
    % Test Z = f(R) for z0 guess
    [Cs,CI]=sort(rStats.r);
    Zmdl = fitlm(rStats.r(CI),rStats.z(CI),'linear','RobustOpts','on');
%     rtest = linspace(0,min(rStats.r),101)';
%     [zpred,zci] = predict(Zmdl,rtest,'Alpha',.01);

    % Get z0 w/ bounds from R = f(z)
    funfun = @(x) getZeroCross(rStats.Rmdl,x);
    rStats.z0 = fminsearch(funfun,Zmdl.Coefficients.Estimate(1));
    funfun = @(x) getLoBoundFun(rStats.Rmdl,x);
    rStats.z0ci(1) = fminsearch(funfun,Zmdl.Coefficients.Estimate(1));
    funfun = @(x) getHiBoundFun(rStats.Rmdl,x);
    rStats.z0ci(2) = fminsearch(funfun,Zmdl.Coefficients.Estimate(1));
    rStats.z0ci = sort(rStats.z0ci);
    
end

% Sub functions to search for zero cross and hi/lo bounds of z0
function a = getZeroCross(Rmdl,x)
%     a = abs(predict(Rmdl,x,'Alpha',.01));
    a = abs(predict(Rmdl,x,'Alpha',.05,'Prediction','observation'));
end

function ci0 = getLoBoundFun(Rmdl,x)
%     [a,b] = predict(Rmdl,x,'Alpha',.01);
    [a,b] = predict(Rmdl,x,'Alpha',.05,'Prediction','observation');
    ci0 = abs(min(b));
end

function ci0 = getHiBoundFun(Rmdl,x)
%     [a,b] = predict(Rmdl,x,'Alpha',.01);
    [a,b] = predict(Rmdl,x,'Alpha',.05,'Prediction','observation');
    ci0 = abs(max(b));
end