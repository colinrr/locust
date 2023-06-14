function Bfit = getTexponent(z,T,idx,z0,zError,r0,sensitivityFlag,T0err,cErr)
  % Power law fitting T decay curve
  % INPUT:
  %     z       = height vector (or time, alternatively)
  %     T       = temperature vector
  %     idx     = optional subset of vector indices
  %     z0      = estimate of virtual origin position (default = 0)
  %     zError  = [lo hi] uncertainty range in z0
  %     r0      = dimensional scale for z vector (default = 1)
  %     T0err   = absolute limits on c coefficient in Kelvin
  %     cErr    = relative limits on c coefficient (fraction of max dT)
if nargin<9 || isempty(cErr)
    cErr = [-0.5 0.5]; % Generously wide defaults
end
if nargin<8 || isempty(T0err)
    T0err = [-40 40]; % Generously wide defaults
end
if nargin<7 || isempty(sensitivityFlag)
    sensitivityFlag = false;
end
if nargin<6 || isempty(r0)
    r0 = 1;
end
if nargin<5 || isempty(zError)
    zError = [NaN NaN];
end
if nargin<4 || isempty(z0)
    z0 = 0;
end
if nargin<3 || isempty(idx)
    idx = 1:length(z);
end

if numel(idx)==2
    idx = idx(1):idx(2);
end

    zError = sort(zError);
    
    Bfit0.zLim = [z(idx(1)) z(idx(end))];
    z = z(idx);
    T = T(idx);
   
    Bfit0.z = z;
    Bfit0.T = T;
    Bfit0.idx = idx;

    Bfit = findB(Bfit0,z0,r0,zError,T0err,cErr);
    
    if sensitivityFlag
        textprogressbar('   Running exponent sensitivity test:  ')
        textprogressbar(0)
        Bfit.z0test = -500:20:min(z); %temp %(2*z0/abs(r0)):.1:(min(z)/abs(r0)))'.*r0;
        Bfit.Btest  = zeros(length(Bfit.z0test),3); % [fitlm fminsrch]
        Bfit.Ctest  = zeros(length(Bfit.z0test),3);
        Bfit.Bnrmse  = zeros(length(Bfit.z0test),1);
        for zi=1:length(Bfit.z0test)
            bb = findB(Bfit0,Bfit.z0test(zi),r0,zError,T0err,cErr);
            % Should also check for values that give theory HERE
            ci = confint(bb.pfit{2});
            Bfit.Btest(zi,:) = [ci(1,2) bb.Bpl(2) ci(2,2)];
            Bfit.Ctest(zi,:) = [ci(1,3) bb.pfit{2}.c ci(2,3)];
            Bfit.Bnrmse(zi,:) = bb.gof{2}.rmse;
%             Bfit.Btest(zi,:) = [bb.B bb.BfSrch1(2)];
%             Bfit.Brms(zi,:)  = [bb.Tmdl.RMSE bb.BfSrch1_rms(2)];
            textprogressbar(zi/length(Bfit.z0test)*100)
        end
    %         [m,i] = min(Brms,[],1);
        
        textprogressbar(' --> Done')
    end
end

function B = findB(B,z0,r0,z0ci,T0err,cErr)
% Estimate B exponent using a few different methods, record results in a
% struct
    T0err = sort(T0err);
    cErr  = sort(cErr);

%     ksdT = linspace(-5,0,151);
    
    normT = @(x) double(x./max(x));
    normZ = @(x,x0,x1) (x-x0)./abs(x1);
    
    T = normT(B.T); %double(B.T./max(B.T));
    z = normZ(B.z,z0,r0);
    B.cLim = [max([T0err(1)./max(B.T) cErr(1)]) min([T0err(2)./max(B.T); cErr(2)])]; % Take minimum error
    B.z0guess   = z0;
    B.z0confInt = z0ci;
    B.r0      = r0;
    
    % Lo, central, hi fit estimates for z0 confidence interval - BEST BY
    % FAR
    [B.pfit{1},B.gof{1}] = fitTpowerLaw(normZ(B.z,z0ci(1),r0),T,B.cLim); 
    [B.pfit{2},B.gof{2}] = fitTpowerLaw(z,T,B.cLim); 
    [B.pfit{3},B.gof{3}] = fitTpowerLaw(normZ(B.z,z0ci(2),r0),T,B.cLim);
    conf = nan(4,3);
    if ~isempty(B.pfit{1}); conf(1:2,:) = confint(B.pfit{1}); end
    if ~isempty(B.pfit{3}); conf(3:4,:) = confint(B.pfit{3}); end
%     if ~isnan(B.pfit{2}.b); b = B.pfit{2}.b; end
%     conf = [confint(B.pfit{1}); confint(B.pfit{3})];
    B.Bpl = [nanmin(conf(1:2,2)) B.pfit{2}.b nanmax(conf(3:4,2))];
    
%% %%%  % Deprecating extra fit methods, they do not perform well %%%%%%%%%
%    % ---- LINEAR FIT MODEL IN LOG-LOG SPACE. LOW CONFIDENCE! ----
%     B.logT = log10(T);
%     B.logZ = log10(z); %(z-z0)./r0);
%     B.dlogT_dz = gradient(B.logT, B.logZ);
% 
%      % REMOVE EXTREME VALUES and get gradient PDF 
%     B.dlogT_dz(B.dlogT_dz>0) = NaN;
%     B.dlogT_dz(B.dlogT_dz<-5) = NaN;
%     B.dTpdf    = ksdensity(B.dlogT_dz,ksdT);
%     [~,ksMaxI]          = max(B.dTpdf);
%     B.dTmode   = ksdT(ksMaxI);  
%     
%     B.Tmdl = fitlm(B.logZ,B.logT,'linear','RobustOpts','on');
%     B.Blm(2) = B.Tmdl.Coefficients.Estimate(2);
%     
%     % B error bounds: LO
%     if ~isnan(z0ci(1))
%         logZ = log10(normZ(B.z,z0ci(1),r0));
%         Tmdl = fitlm(logZ,B.logT,'linear','RobustOpts','on');
%         Bci = coefCI(Tmdl,.01);
%         B.Blm(1) = min(Bci(2,:));
%     else
%         B.Blm(1) = NaN;
%     end
%     
%     % B error bounds: HI
%     if ~isnan(z0ci(2))
%         logZ = log10(normZ(B.z,z0ci(2),r0));
%         Tmdl = fitlm(logZ,B.logT,'linear','RobustOpts','on');
%         Bci = coefCI(Tmdl,.01);
%         B.Blm(3) = max(Bci(2,:));
%     else
%         B.Blm(3) = NaN;
%     end
%    
%    % ------------------------------------------------------------
%
%     % ---- FMINSEARCH METHOD: MINIMIZE LINEAR ERROR FOR NORMALIZED T vs z^B ----
%     % - ALSO LOW CONFIDENCE! - better than lin-log fit though
%     normfun = @(x) (x-min(x))./range(x); % Normalization
%     
%     % Fminsearch method, find B, fixed z0 w/ lo/hi bounds
% %     funner = @(x) rms(T-((B.z-z0).^x./max((B.z-z0).^x)));
%     funner = @(x) rms(normfun(B.T)-normfun((B.z-z0ci(1)).^x));
%     B.BfSrch1(1) = fminsearch(funner,B.Blm(2));
%     B.BfSrch1_rms(1) = funner(B.BfSrch1(1));
%     
%     funner = @(x) rms(normfun(B.T)-normfun((B.z-z0).^x));
%     B.BfSrch1(2) = fminsearch(funner,B.Blm(2));
%     B.BfSrch1_rms(2) = funner(B.BfSrch1(2));
%     
%     funner = @(x) rms(normfun(B.T)-normfun((B.z-z0ci(2)).^x));
%     B.BfSrch1(3) = fminsearch(funner,B.Blm(2));
%     B.BfSrch1_rms(3) = funner(B.BfSrch1(3));
%     
%     [B.BfSrch1,bi] = sort(B.BfSrch1);
%     B.BfSrch1_rms = B.BfSrch1_rms(bi);
%   
%%%%%%%  % ------------------------------------------------------------

end

% Power law fit function
function [fitresult, gof] = fitTpowerLaw(z,T,cLim)
 % C lim = optional [low high] bounds on c coefficient 

    if nargin<3 || isempty(cLim)
        cLim = [-0.5 0.5];
    end

    if z(1)<=0
        warning('Virtual origin gives negative z values!')
%         fitresult   = cfit(fittype('a*x^b+c'),nan,nan,nan);
        fitresult   = [];
        gof.sse     = nan;
        gof.rsquare = nan;
        gof.dfe     = nan;
        gof.adjrsquare = nan;
        gof.rmse       = nan;
    else
        
        % Quick and dirty to explore weighting the earliest points
        weights = ones(size(T));  % Default unity weights
%         weights(1:round(length(T)*0.2)) = 5; % test weight

        [xData, yData, weights] = prepareCurveData( z, T, weights );

        % Set up fittype and options.
        ft = fittype( 'power2' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
%         opts.Robust = 'Bisquare';
%         opts.Robust = 'LAR';
        opts.StartPoint = [1./(z(1).^(-5/3)) -5/3 0];
        opts.Lower = [0 -10 min(cLim)];
        opts.Upper = [Inf 0 max(cLim)];
        opts.Weights = weights;

        [fitresult, gof] = fit( xData, yData, ft, opts );
    end
    
end