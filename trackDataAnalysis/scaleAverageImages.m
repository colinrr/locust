
% Get estimates of  r0,z0 from avg Img and track bounds
clear all

dataDir   = '~/Kahuna/data/sabancaya_5_2018/';
% dataDir = 'C:\Users\crowell\Kahuna\data\sabancaya_5_2018';

maskFile = fullfile(dataDir,'pulseTrack_analysis/allMasks_3events.mat');  % Plume masks only

% avgImgFile = fullfile(dataDir,'pulseTrack_analysis/eventAveragedImages_2022-07-12.mat');
% avgImgFile = fullfile(dataDir,'pulseTrack_analysis/eventAveragedImages_2022-07-16.mat');
% avgImgFile = fullfile(dataDir,'pulseTrack_analysis/eventAveragedImages_2022-11-02.mat');
% avgIdc     = 1:5; % Which average images to use?
% 
% trackPolyFile = fullfile(dataDir,'pulseTrack_analysis/tracks_imgs_polys_n26_2022-07-04.mat');
% tkIdc = [4	84	23	422];
% 
% nS = 10;
% 
% % Indices [rIdx1 rIdx2] for drdz,z0 estimation - [gaussian mask bwdist] 
% % rIdxList = {[1 162] [1 167] [1 153];
% %             [1:103 168:202] [1:102 181:202] [1 169];
% %             [1 234] [1 234] [1 234];
% %             [1 223] [1 170] [1 166];
% %             [1 157] [1 244] [1 177]};
%         
% % NEW Indices [rIdx1 rIdx2] for drdz,z0 estimation - [gaussian mask VzGauss] 
% rIdxList = {[1 362] [1 167] [1 147];
%     [1:103 168:202] [1:102 181:202] [1 262];
% %     [1:103] [1:102] [1 262];
%     [1 234] [1 130] [1 130];
%     [1 258] [1 258] [1 250];
%     [1 267] [1 261] [1 267]};
% 
% rAvgIdx = {[1:147],[1:102 181:202],[1:130],[1:250], [1:261]};
% % rAvgIdx = {[1:147],[1:102],[1:130],[1:250], [1:261]};
% 
% TIdxList = {[],[],[],[1 292],[1 298]};

% ##### NEW AVG IMG FILE PARAMS - REORDERED
avgImgFile = fullfile(dataDir,'pulseTrack_analysis/eventAveragedImages_2022-11-02.mat');
avgIdc     = 1:6; % Which average images to use?

% trackPolyFile = fullfile(dataDir,'pulseTrack_analysis/tracks_imgs_polys_n26_2022-07-04.mat');
trackPolyFile = fullfile(dataDir,'pulseTrack_analysis/tracks_gauss_polys_n26_2022-12-05.mat');
tkIdc = [4	84	23	422];

nS = 10;

% Indices [rIdx1 rIdx2] for drdz,z0 estimation - [gaussian mask bwdist] 
% rIdxList = {[1 162] [1 167] [1 153];
%             [1:103 168:202] [1:102 181:202] [1 169];
%             [1 234] [1 234] [1 234];
%             [1 223] [1 170] [1 166];
%             [1 157] [1 244] [1 177]};
        
% NEW Indices [rIdx1 rIdx2] for drdz,z0 estimation - [gaussian mask VzGauss] 
rIdxList = {[1 234] [1 130] [1 130]; % 25A
    [1:103 168:202] [1:102 181:202] [1 262]; % 24A
%     [1:103] [1:102] [1 262];
    [1 362] [1 167] [1 147]; % 25B main pulses
    [1 258] [1 258] [1 250]; % 25B full history
    [1 145] [1 153] [1 171];              % 25B first pulse
    [1 267] [1 261] [1 267]}; % 25B decay period

rAvgIdx = {[1:130],[1:102 181:202],[1:147],[1:250],[1 152], [1:261]};
% rAvgIdx = {[1:147],[1:102],[1:130],[1:250], [1:261]};

TIdxList = {[],[],[],[1 292],[],[1 298]};

% #############################

sensitivityFlag = false;
dT0Limits       = [-10 10]; % Absolute c error limits about 0 (K)
cLimits         = [-0.2 0.2]; % Relative c error limits about 0

outFile = fullfile(dataDir,'pulseTrack_analysis/AverageImageScalingResults_V3.mat');
saveOutput = false;
%%
% load(trackPolyFile)r
% load(maskFile)
load(avgImgFile)

%%

% Output data summary
z0      = zeros(5,3,4);
r0      = zeros(5,3,4);
drdz    = zeros(5,3,4);
B       = zeros(5,3,2);

for ii = avgIdc
    ridx = rIdxList(ii,:);
    Tidx = TIdxList{ii};
    IS = Tavg(ii);
    
    IS.rGstats = getRparams(IS.z,IS.gaussCoeffs(:,3),ridx{1});  % radius,drdz,z0 from gaussian profile fits
    IS.rMstats = getRparams(IS.z,IS.r_mask,ridx{2});  % radius,drdz,z0 from mask widths
%     IS.rBWstats = getRparams(IS.z,double(IS.r_bw),ridx{3});  % radius,drdz,z0 from mask bwdist calculation
    IS.rVstats = getRparams(IS.z,IS.VzGaussCoeffs(:,3),ridx{3}); % radius,drdz,z0 from Vz gaussian profile fits
    
    % Averaging the radii params
%     weights = [IS.rGstats.Rmdl.Rsquared.Adjusted IS.rMstats.Rmdl.Rsquared.Adjusted IS.rVstats.Rmdl.Rsquared.Adjusted].^2;
    weights = [1 1 1];
%     if ii==1 % Checking that these don't matter because they are scalar multiples
%         rScales = [sqrt(2) sqrt(2) 0.71]; % Patrick 2007b summary of conversions to top-hat radius for plumes
%     else
%         rScales = [sqrt(2) sqrt(2) 0.9]; % Patrick 2007b summary of conversions to top-hat radius for thermals
%     end
%     rAll = [IS.gaussCoeffs(:,3) IS.VzGaussCoeffs(:,3) IS.r_mask] .* rScales; 
    rAll = [IS.gaussCoeffs(:,3) IS.VzGaussCoeffs(:,3) IS.r_mask];
    IS.rAll = sum(rAll.*weights,2)./sum(weights);
%     IS.rAll = mean([IS.gaussCoeffs(:,3) IS.VzGaussCoeffs(:,3) IS.r_mask],2);
    IS.rAstats = getRparams(IS.z,IS.rAll,rAvgIdx{ii});
    
%     Tz = prctile(IS.Timg,95,2);
    Tz = IS.Tz; 
    IS.Bfit_rG = getTexponent(IS.z,Tz,Tidx,...
        IS.rGstats.z0,IS.rGstats.z0ci,IS.rGstats.r0,sensitivityFlag,dT0Limits,cLimits);  % T exponent fit using gaussian radius fitting

    IS.Bfit_rM = getTexponent(IS.z,Tz,Tidx,...
        IS.rMstats.z0,IS.rMstats.z0ci,IS.rMstats.r0,sensitivityFlag,dT0Limits,cLimits);  % T exponent fit using gaussian radius fitting

%     IS.Bfit_rBW = getTexponent(IS.z,IS.Tz,Tidx,...
%         IS.rBWstats.z0,IS.rBWstats.z0ci,IS.rBWstats.r0,sensitivityFlag);  % T exponent fit using gaussian radius fitting

    IS.Bfit_rV = getTexponent(IS.z,Tz,Tidx,...
        IS.rVstats.z0,IS.rVstats.z0ci,IS.rVstats.r0,sensitivityFlag,dT0Limits,cLimits);  % T exponent fit using gaussian radius fitting

    % ###### METHOD OPTIONS: GETTING A COMBINED FIT ESTIMATE ######
    % (1)-- Use the FULL range of z0 estimates from different methods
%     z0vec = [IS.rGstats.z0 IS.rMstats.z0];% IS.rBWstats.z0];
%     z0vec = [IS.rGstats.z0 IS.rMstats.z0 IS.rVstats.z0];
%     z0mu = mean(z0vec);
%     z0ci = [min(z0vec) max(z0vec)];
%     r0mu = mean([IS.rGstats.r0 IS.rMstats.r0]);% IS.rBWstats.r0]);

    % (2)-- Fitting to an average radius
    z0mu = IS.rAstats.z0;
    z0ci = IS.rAstats.z0ci;
    r0mu = IS.rAstats.r0;
    
    % (3)-- Using the mean z0 AND mean confidence interval
%     z0rng = [IS.rGstats.z0ci(1) IS.rMstats.z0ci(1) IS.rVstats.z0ci(1); % Mean low bound
%              IS.rGstats.z0 IS.rMstats.z0 IS.rVstats.z0; % Mean central estimate
%              IS.rGstats.z0ci(2) IS.rMstats.z0ci(2) IS.rVstats.z0ci(2)]; % Mean upper bound
%     z0mu  = mean(z0rng(2,:));
%     z0ci = mean(z0rng([1 3],:),2);
%     r0mu = mean([IS.rGstats.r0 IS.rMstats.r0 IS.rVstats.r0]);
    
    % -- Get the fit
    IS.Bfit = getTexponent(IS.z,Tz,Tidx,...
        z0mu,z0ci,r0mu,sensitivityFlag,dT0Limits,cLimits);

    TavgS(ii) = IS;
    
    % Output table
    z0(ii,:,1) = sort([IS.rGstats.z0ci IS.rGstats.z0]);
    z0(ii,:,2) = sort([IS.rVstats.z0ci IS.rVstats.z0]);
    z0(ii,:,3) = sort([IS.rMstats.z0ci IS.rMstats.z0]);
    z0(ii,:,4) = sort([IS.rAstats.z0ci IS.rAstats.z0]);
    r0(ii,:,1) = sort([IS.rGstats.r0ci IS.rGstats.r0]);
    r0(ii,:,2) = sort([IS.rVstats.r0ci IS.rVstats.r0]);
    r0(ii,:,3) = sort([IS.rMstats.r0ci IS.rMstats.r0]);
    r0(ii,:,4) = sort([IS.rAstats.r0ci IS.rAstats.r0]);
    drdz(ii,:,1) = sort([IS.rGstats.drdzCI IS.rGstats.drdz]);
    drdz(ii,:,2) = sort([IS.rVstats.drdzCI IS.rVstats.drdz]);
    drdz(ii,:,3) = sort([IS.rMstats.drdzCI IS.rMstats.drdz]);
    drdz(ii,:,4) = sort([IS.rAstats.drdzCI IS.rAstats.drdz]);
    B(ii,:,1) = IS.Bfit_rG.Bpl;
    B(ii,:,2) = IS.Bfit.Bpl;
end
%     rMaskStats = getRparams(Img.z,Img.r_mask,IdcM);
%     rBWstats   = getRparams(Img.z,Img.r_bw,IdcBW);

% Output table
varNames = {'z0_Tgauss','z0_Vgauss','z0_mask','z0_combo','r0_Tgauss','r0_Vgauss','r0_mask','r0_combo',...
    'drdz_Tgauss','drdz_Vgauss','drdz_mask','drdz_combo','B_Tgauss','B_combo'};
AvgImgT = table(z0(:,:,1),z0(:,:,2),z0(:,:,3),z0(:,:,4),r0(:,:,1),r0(:,:,2),r0(:,:,3),r0(:,:,4),...
    drdz(:,:,1),drdz(:,:,2),drdz(:,:,4),drdz(:,:,4),B(:,:,1),B(:,:,2),...
    'RowNames',strrep({TavgS.description},' ','_'),'VariableNames',varNames);
%% Plot averaged image results
dx = 0.01;
dx2 = 0.05;
dy = 0.1;
ppads = [0.05 0.02 0.12 0.08];
nr = 1; %length(avgIdc);
nc = 4;
lw = 1.2;
xSz = [2 1 1 3];
fs = 14;
 
% Img, r(z) + z0
for ii = 1:length(avgIdc)
    
    IS = TavgS(avgIdc(ii)); % TavgS(avgIdc(ii)).Istats;
    
    %     Bfits = [IS.Bfit_rG IS.Bfit_rM IS.Bfit_rBW IS.Bfit];
    Bfits = [IS.Bfit_rG IS.Bfit_rM IS.Bfit_rV IS.Bfit];
    Tfit = cell(1,length(Bfits));
    BAll = [];
    for jj = 1:length(Bfits)
%         BAll = [BAll Bfits(jj).Bpl(2) Bfits(jj).BfSrch1(2) nanmean(Bfits(jj).dlogT_dz)];
        [Tfit{jj},~] = getDimensionalTcurve(Bfits(jj),Bfits(jj).z);
%         Tfit{jj} = Td{2};
%         Tfit{jj} = Td{2}.*max(Bfits(jj).T);
    end
%     BmeanAll = mean(BAll);

    figure('position',[100 200 1400 800],'name',IS.description)
    pi = nc*(ii-1) + 1;
    
    % IMAGE
    ax(1) = tightSubplot(nr,nc,1,dx,dy,ppads,xSz); 
    plotThermVelocities(IS.x,IS.z,[],[],'thermal',IS.Timg,'mask',IS.mask(IS.maskI,:),...
        'axis',gca,'Trange',[min(IS.Timg(:)) max(IS.Timg(:))]);
    yl = ylim;
    
    % ####### r(z) and z0 estimates ######
    ax(2) = tightSubplot(nr,nc,2,dx,dy,ppads,xSz); 
    co = get(gca,'ColorOrder');
    co = co([1 2 4:7],:);
    hold on
    plot(IS.gaussCoeffs(:,3),IS.z,':','Color',rgba2rgb(co(ii,:),0.5),'LineWidth',lw)
    plotZ0estimate(IS.rGstats,ax(2),co(1,:),lw)
    lp(1) = plot(nan,nan,'Color',co(1,:),'LineWidth',lw);
    
    
    plot(IS.r_mask,IS.z,':','Color',rgba2rgb(co(2,:),0.5),'LineWidth',lw)
    plotZ0estimate(IS.rMstats,ax(2),co(2,:),lw)
    lp(2) = plot(nan,nan,'Color',co(2,:),'LineWidth',lw);
    
%     plot(IS.r_bw,IS.z,':','Color',rgba2rgb(co(3,:),0.5),'LineWidth',lw)
%     plotZ0estimate(IS.rBWstats,ax(2),co(3,:),lw)
%     lp(3) = plot(nan,nan,'Color',co(3,:),'LineWidth',lw);

    plot(IS.VzGaussCoeffs(:,3),IS.z,':','Color',rgba2rgb(co(3,:),0.5),'LineWidth',lw)
    plotZ0estimate(IS.rVstats,ax(2),co(3,:),lw)
    lp(3) = plot(nan,nan,'Color',co(3,:),'LineWidth',lw);

    plot(IS.rAll,IS.z,':k','LineWidth',lw)
    plotZ0estimate(IS.rAstats,ax(2),[0 0 0],lw)
    lp(4) = plot(nan,nan,'Color','k','LineWidth',lw);
    
    lp(5) = errorbar(0,IS.Bfit.z0guess,IS.Bfit.z0guess-IS.Bfit.z0confInt(1),IS.Bfit.z0confInt(2)-IS.Bfit.z0guess,...
        'o','Color',co(4,:),'MarkerSize',10,'LineWidth',1.8);
    
%     lg = legend(lp,{'Gaussian','Mask width','BW dist'});
    lg = legend(lp,{'r(T_{Gauss})','r(mask)','r(w_{Gauss})','r(Averaged)','Combined z0'});
    axis equal
    xlim([-10 yl(2)/2])
    grid on
    xlabel('r (m)')
    
    % T(z), with fit selection and curve confidence bounds
    ax(3) = tightSubplot(nr,nc,3,dx,dy,ppads,xSz); % T(z)
    plot(IS.Tz,IS.z,':','Color',[0.5 0.5 0.5],'LineWidth',lw)
    hold on
    plot(IS.Bfit_rG.T,IS.Bfit_rG.z,'.','Color','k','LineWidth',lw)
    % T fit confidence intervals are essentially identical for all choices
    % of z0 - just the coeffs change
    for jj = 1:length(Bfits)
        plot(Tfit{jj}(:,2),IS.Bfit_rG.z,'Color',co(jj,:),'LineWidth',lw)
        plot(Tfit{jj}(:,[1 3]),IS.Bfit_rG.z,'--','Color',co(jj,:),'LineWidth',lw)
    end
    grid on
    xlabel('\Delta T')
    
    linkaxes(ax(1:3),'y')
    axis(ax(2),'tight')
    
    % Power law fits
    ax(4) = tightSubplot(nr*2,nc,4,dx2,dy,ppads,xSz); % dr/dz estimates
    ax(5) = tightSubplot(nr*2,nc,8,dx2,dy,ppads,xSz); % B estimates
    

%     plotBestimate(Bfits,ax(4:5),{'Gauss Fit','Mask width','BW dist','Combined'},co([1 2 3 4],:))
    plotBestimate(Bfits,ax(4:5),{'T Gauss','Mask width','Vz Gauss','Combined'},co([1 2 3 4],:))
    xl4 = xlim(ax(4));
%     plot(ax(4),xl4,BmeanAll.*[1 1],'Color',[0.5 0.5 0.5],'LineWidth',2)
    ll = get(ax(4),'Legend');
    ll.String{jj+1} = 'Mean Estimate';

    set(ax,'FontSize',fs)
end

%%
% Choose output subset and rename


if saveOutput
    save(outFile,'TavgS','AvgImgT')
end
%% QUICK PLOTS TO QC CURVE FITS

% r(z) curves vs index
% for ii = 1:5
%     IS = Tavg(ii);
%     figure; 
%     subplot(1,4,1); 
%     imagesc(Tavg(ii).VzMasked); set(gca,'YDir','normal');colorbar;colormap(plasmagrey)
%     cax = caxis;
%     subplot(1,4,2)
%     imagesc(Tavg(ii).VzSynth); set(gca,'YDir','normal');colorbar;colormap(plasmagrey)
%     caxis(cax)
%     subplot(1,4,3)
%     zi = 1:length(IS.z);
%     plot(IS.gaussCoeffs(:,3),zi);hold on; plot(IS.VzGaussCoeffs(:,3),zi);
%     plot(IS.r_mask,zi,'k')
%     legend({'r_{Tgauss}','r_{mask}','r_{Wgauss}'})
%     subplot(1,4,4)
%     plot(IS.Grmse./max(IS.Timg,[],2),zi);hold on; plot(IS.VzGrmse./max(IS.VzMasked,[],2),zi);
%     legend('T Gaussian RMSE','W Gaussian RMSE')
% end
%%
% trackSet = [1 9 20]; %[1 5 9 19 20 24 26];

% for ti = trackSet
%     x    = tk(ti).x - tk(ti).x0;
%     z    = tk(ti).z;
%     mask = tk(ti).trackBounds.mask;
% 
%     dx = mean(diff(z));
% %     pixDist = bwdist(~mask);
% 
%     % x_center = zeros(length(z));  
% 
%     % Get mask skeleton
%     % BW3 = bwmorph(mask,'skel',Inf);
% %     BW = bwskel(mask,'MinBranchLength',ceil(max(pixDist(:)))); % In theory should get 1 single branch
% %     idx = find(BW);
% % 
% %     % Get pixel indices and sort by distance from vent
% %     [ii,jj] = ind2sub(size(BW),idx);
% %     s = sqrt(x(jj).^2 + z(ii).^2);
% %     [s,si] = sort(s);
% %     ii = ii(si);
% %     jj = jj(si);
% %     idx = idx(si);
% %     rpix = pixDist(idx);
% %     r = rpix.*dx;
%     rnpx = (tk(ti).npx./pi).^.5.*dx; % Cluster-area R
%     
%     tI = 1:1:length(tk(ti).clustZ);
%     cx = zeros(length(tI),1);
%     cz = cx;
%     for kk=1:length(tI)
%         [cii,cjj] = ind2sub(size(mask),tk(ti).clustIdx{tI(kk)});
%         cx(kk) = mean(x(cjj));
%         cz(kk) = mean(z(cii));
%     end
%     s = sqrt(cx.^2 + cz.^2);
%     dsdt = gradient(s,tk(ti).t); % Velocity
%     nS  = 30; %rnpx(1)/mean(dsdt)/mean(diff(tk(ti).t))
%     cxS = smooth(cx,nS);
%     czS = smooth(cz,nS);
%     
%     dzdx = gradient(czS,cxS);
% %     angle = atan2d(czS,cxS); % Angle to horizontal of track axis
%     % Coeffs for perpendicular lines
%     dzdx2 = -1./(dzdx);
%     b     = czS - dzdx2.*cxS;
%     xp = x([min(tk(ti).trackBounds.X) max(tk(ti).trackBounds.X)])';
%     zp = dzdx2.*xp + b;
%     tbS.X = smooth(x(tk(ti).trackBounds.X),nS/2);
%     tbS.Z = smooth(z(tk(ti).trackBounds.Y),nS/2);
%     
%     PXint = zeros(size(zp));
%     PZint = zeros(size(zp));
%     for xi=1:size(zp,1)
% %         P = InterX([xp; zp(xi,:)],[x(tk(ti).trackBounds.X) , z(tk(ti).trackBounds.Y)]');
%         P = InterX([xp; zp(xi,:)],[tbS.X , tbS.Z]');
%         PXint(xi,:) = P(1,[1 end]);
%         PZint(xi,:) = P(2,[1 end]);
%     end
%     r = sqrt( ( cxS - PXint ).^2 + ( czS - PZint ).^2 );
%     rmu = mean(r,2);
%     %%
%     figure('position',[50 50 1000 800])
% %     subplot(2,2,1)
% %     imagesc(x,z,pixDist)
% %     set(gca,'YDir','normal')
% %     colorbar
% %     hold on
% %     plot(x(jj),z(ii))
% %     plot(x(tk(ti).trackBounds.X),z(tk(ti).trackBounds.Y))
% %     colormap(gray)
%     
%     subplot(2,2,2)
%     plotTrackOutlines(tk(ti),1)
%     hold on
% 
%     set(gca,'ColorOrderIndex',1)
%     plot(cx+tk(ti).x0,cz,'.')
%     plot(cxS+tk(ti).x0,czS,'LineWidth',1.5)
% 
%     subplot(2,2,3)
%     hold on
%     plot(czS,r,'.')
%     plot(czS,rmu,'k')
%     plot(tk(ti).clustZ,rnpx,'.')
%     plot(tk(ti).Istats.z,tk(ti).Istats.gaussCoeffs(:,3),'.')
% 
%     subplot(2,2,4)
%     plot(s,rnpx,'.')
% end