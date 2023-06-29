function [D,atmo,fitStats] = fitAtmoProfile(D,atmo,tI,imgROI,TErrThresh,zUncThresh,TmaskThresh)
% D = fitAtmoProfile(D,atmo,tI,zI,TErrorThresh,zUncThresh)
% Fit atmospheric profiles to thermal imagery, estimate Temperature DC shift
% for a robust thermal data threshold.
% IN:
%   D           = data cube struct or path to data cube struct
%   atmo        = atmospheric profile cell of structs from "getAtmoProfile.m"
%   tI          = vector frame subscripts to subset D.T [Default ALL]
%   imgROI      = [I1 I2 J1 J2] Subscripts within image frame, excludes 
%                 pixels outside this region. [Default ALL]
%   REMOVED: zI          = vector of height subscripts to subset D.T [Default ALL]
%
%   The following uncertainty thresholds are used to crop out pixels that
%   are likely to have significant error in height estimate:
%
%   TErrThresh = Temperature error used to calculate height error threshold
%                   (error due to 3D SHAPE OF PLUME)
%                 Calculated according to:  [DEFAULT = 1.5 K]
%                  zError = TErrorThresh./(Mean profile lapse rate)
%                   
%   zUncThresh  = Z uncertainty threshold for height errors resulting from
%               unknown PLUME AXIAL POSITION into/out of the viewing plane.
%               (Distinct from Z error above). Typically threshold is 
%               similar to or larger than zError. [DEFUALT = zError*1.5]
%
%   TmaskThresh = Temperature threshold to use for extra mask clippingsubsetI.
%               This helps removes background pixels for things likes holes
%               in the plume. Pick a value from near the plume edge in
%               colder upper parts of the plume, where the temperature
%               clearly starts to fall off from plume towards background.
%
% OUT:
%   D       = data cube struct with updated atmospheric profile paramters
%               included.
%   atmo    = updated atmo struct (or cell of structs for multiple
%               profiles). If atmo input is entered as a path, will update
%               the saved file with calculated fit and error parameters.
%
% C Rowell Jun 2020
fprintf('\n========= FIT Atmospheric Profile(s) to Thermal Data =========\n')

if nargin<7
    TmaskThresh = [];
end
if nargin<6
    zUncThresh = [];
end
if nargin<5
    TErrThresh = [];
end
if nargin<4
    imgROI = [];
end
if nargin<3
    tI = [];
end
%% data prep
if ischar(D)
    if exist(D,'file')
        disp('Loading data cube...')
        load(D)
    end
end
if ischar(atmo)
    if exist(atmo,'file')
        atmoFile = atmo;
        load(atmo)
    end
else
    atmoFile = [];
end
if isempty(tI)
    tI = 1:size(D.T,3);
end
if isempty(imgROI)
    zI = 1:size(D.T,1);
    xI = 1:size(D.T,2);
else
    assert(length(imgROI)==4,'imgROI format: [i1 i2 j1 j2]')
    zI = imgROI(1):imgROI(2);
    xI = imgROI(3):imgROI(4);
end

clippedMaskAll = D.mask;
if ~isempty(TmaskThresh)
    clippedMaskAll(D.T<TmaskThresh) = false;
end
for ll=1:size(clippedMaskAll,3)
    clippedMaskAll(:,:,ll) = bwdist(~clippedMaskAll(:,:,ll)) > 20;
end
clippedMask = clippedMaskAll(zI,xI,tI);

varioCamTres = 0.05; % VarioCam temperature resolution
TresStep     = 4; % Multiplier of resolution to get PDF temperature step
dT           = varioCamTres*TresStep;
%% SET UP LOOP

    % Raw temperature data, plume edges cut out
    Z = D.z(zI);
    X = D.x(xI);
    tt = D.t(tI);
    
    % z Error estimation
    disp('Estimating height error...')
    [zErr,zErrMax] = zErrorEstimation(X,Z+(D.geom.Ztarg-D.geom.Z0),clippedMask,D.geom,'cylindrical');
    zErrUncertainty = squeeze(zErrMax(:,3,:)-zErrMax(:,1,:));

    rmsFig = figure;
    ncols = length(atmo)+1;
    rmsAx = subplot(1,ncols,ncols);
    
    disp('Calculating profile fits...')
    for kk = 1:length(atmo)
        atmo{kk}.subsetI = zI;
        atmo{kk}.subsetK = tI;
        
        rawT = (D.T - atmo{kk}.Tinterp).*clippedMaskAll;
        rawT(rawT==0) = NaN;
        % Remove profile and apply NaNs
        clippedT = rawT(zI,xI,tI);
        % Temperature vector for PDFs
        Tlims = prctile(rawT(:),[1 99]);
        ksT = round(nanmin(rawT(:))./(dT))*dT:dT:round(nanmax(rawT(:))./(dT))*dT;
        
        dTmin_crop = squeeze(nanmin(clippedT,[],2)); % Apparent minima of subset
        % Probability distrubution of minima over time
        pdmin_crop = ksdensity(dTmin_crop(:),ksT);

        % ------ ERROR THRESHLOLDS -------
        if isempty(TErrThresh)
            atmo{kk}.TErrThresh = 1.5;
            % Approximate height difference for 1.5 K uncertainty
        else
            atmo{kk}.TErrThresh = TErrThresh;
        end
        atmo{kk}.zErrThresh = abs([atmo{kk}.TErrThresh]/mean(gradient(atmo{kk}.Tinterp,mean(diff(Z)))));
        fprintf('Calculated z error threshold: +/- %.2f m\n',atmo{kk}.zErrThresh)
        
        if isempty(zUncThresh)
            % 2 K margin for positional uncertainty
            atmo{kk}.zUncThresh = atmo{kk}.zErrThresh*1.5; %abs([2]/mean(gradient(atmo{kk}.Tinterp,mean(diff(D.z)))));
        else
            atmo{kk}.zUncThresh = zUncThresh;
        end
        
        % Get error threshold flags
        zErrFlags = zErr>atmo{kk}.zErrThresh;
        zUncFlags = zErrUncertainty>atmo{kk}.zUncThresh;
        zUncFlags = repmat(permute(zUncFlags,[1 3 2]),[1 size(clippedT,2) 1]);
        
        % Percentages of pixels over thresholds at each height
        zErrN   = sum(zErrFlags,[2 3],'omitnan');
        zErrPct = 1 - zErrN./sum(clippedMask,[2 3],'omitnan');
        zUncN   = sum(zUncFlags.*clippedMask,[2 3],'omitnan');
        zUncPct = 1 - zUncN./sum(clippedMask,[2 3],'omitnan');
        zAllPct = sum(~zErrFlags.*~zUncFlags.*clippedMask,[2 3],'omitnan')./sum(clippedMask,[2 3],'omitnan');
        zNpixFilt = sum(~zErrFlags.*~zUncFlags.*clippedMask,[2 3],'omitnan');
        zFiltFlag  = zNpixFilt>100;
        ErrNaN = zAllPct==0;
        
        % Apply error threshold flags
        Tcut = clippedT;
        Tcut(zErrFlags) = NaN;
        Tcut(zUncFlags) = NaN;
        % ------ -------------- -------
        
        % Distributions over height, after cutting times, heights, and zErrors
        [pd_vs_z_zErrCut_tcut,~,pd_zErrCut_tcut,stats_Tcut] = ksdensityND(Tcut,ksT,[2 3]);
        dTmin_zErrCut = squeeze(nanmin(Tcut,[],2)); % Apparent minima of subset
        

        % Get lower half-max of modal peak
        [pdMax,pdmaxI] = max(pd_zErrCut_tcut);
        % Half maxima
        pdHalfMax = pdMax/2;
        [~,hi] = closest(pdHalfMax,pd_zErrCut_tcut(1:pdmaxI),1);
        
        atmo{kk}.Tmode = ksT(pdmaxI); % TEMPERATURE MODE OF ALL FILTERED DATA
        atmo{kk}.T_halfMax = min(ksT(hi));

        % Fit and residuals
        Tout = Tcut - atmo{kk}.Tmode;
        % mean Sum of squares
        atmo{kk}.Tssq = sum((Tcut(~ErrNaN,:,:)-atmo{kk}.Tmode).^2,[2 3],'omitnan')./sum(clippedMask(~ErrNaN,:,:),[2 3],'omitnan');
        atmo{kk}.zFit = Z(~ErrNaN);
        atmo{kk}.Trms = rms(Tcut(~isnan(Tcut(:)))-atmo{kk}.Tmode);
        atmo{kk}.modeTrms = rms(stats_Tcut.Mode(zFiltFlag)-atmo{kk}.Tmode);
        %% Uncomment this part for extra checks?
        dTmin_all  = squeeze(nanmin(rawT,[],2));

        % Distributions for all height over time
        [pd_vs_z,~,pdall,stats_vs_z_all] = ksdensityND(rawT,ksT,[2 3]);
        rawTrms = rms(stats_vs_z_all.Mode-atmo{kk}.Tmode);
        % Probability distrubution of minima over time
        pdmin = ksdensity(dTmin_all(:),ksT);
        pdmin_zErrCut = ksdensity(dTmin_zErrCut(:),ksT);
        pdmin_all = ksdensity(dTmin_all(:),ksT);
    
         pd_vs_z_crop = ksdensityND(clippedT,ksT,[2 3]);
       
        %% PLOTS

        % -------------- PDF's, fit, error vs height ------------------
        % Plot CROPPED time distrubtion vs height with zError cut
        figure('position',[20 50 1600 600],'name',sprintf(atmo{kk}.quickName))
        co = get(gca,'ColorOrder');

        ax(1) = subplot(1,4,1);
        pcolor(ksT,D.z,(pd_vs_z)')
%         pcolor(ksT,D.z,(pd_vs_z./max(pd_vs_z,[],1))')
        shading flat
        colormap(gray(150))
        caxis([0 1])
        hold on
        pm = plot(stats_Tcut.Mode,Z,'.','Color',co(2,:));
        xlabel('T_{camera} - T_{profile} [K]')
        ylabel('z [m above vent]')
        xlim(prctile(rawT(:),[2 95]))
        legend(pm,{'Filtered T_{Mode}(z)'})
        title('PDF: T(z)_{camera} - T(z)_{profile}, All Pixels')
        set(gca,'FontSize',12)

        ax(2) = subplot(1,4,2);
        pcolor(ksT,Z,(pd_vs_z_crop)')
%         pcolor(ksT,Z,(pd_vs_z_crop./max(pd_vs_z_crop,[],1))')
        shading flat
        colormap(gray(150))
        caxis([0 1])
        hold on
        pm = plot(stats_Tcut.Mode,Z,'.','Color',co(2,:));
        xlabel('T_{camera} - T_{profile} [K]')
        ylabel('z [m above vent]')
        xlim(prctile(rawT(:),[2 95]))
        title('PDF: T(z)_{camera} - T(z)_{profile}, Filtered: t, z')
        legend(pm,{'Filtered T_{Mode}(z)'})
        set(gca,'FontSize',12)

        ax(3) = subplot(1,4,3);
        pcolor(ksT,Z,pd_vs_z_zErrCut_tcut')
%         pcolor(ksT,Z,(pd_vs_z_zErrCut_tcut./max(pd_vs_z_zErrCut_tcut,[],1))')
        shading flat
        colormap(gray(150))
        caxis([0 1])
        hold on
        pk(1) = plot(stats_Tcut.Mode,Z,'.','Color',co(2,:));
        pk(2) = plot(atmo{kk}.Tmode*[1 1],[min(Z) max(Z)],'Color',co(1,:),'LineWidth',1.5);
        pk(3) = plot(atmo{kk}.T_halfMax*[1 1],[min(Z) max(Z)],'--','Color',co(1,:),'LineWidth',1.5);
        xlabel('T_{camera} - T_{profile} [K]')
        ylabel('z [m above vent]')
        xlim(prctile(rawT(:),[2 95]))
        title('PDF: T(z)_{camera} - T(z)_{profile}, Filtered: t, z, zError')
        legend(pk,{'Filtered T_{Mode}(z)','Final T_{Mode}','T_{1/2 Max}'})
        set(gca,'FontSize',12)
        
        ax(4) = subplot(1,4,4);
        plot([zErrPct zUncPct zAllPct],Z,'LineWidth',1.7)
        xlabel('Fraction of pixels below cutoff')
        ylabel('z [m above target]')
        legend('Z Error (plume thickness)','Z Error (plume position)','All filters')
        title('Height uncertainty filter')
        linkaxes(ax,'y')
        set(gca,'FontSize',12)
        
        % --------- Plot distribution functions (total and minima) --------
        figure('name',atmo{kk}.quickName,'position',[100 100 1000 400])
        plot(ksT,pdall,'LineWidth',1.5)
        hold on
        plot(ksT,pd_zErrCut_tcut,'LineWidth',1.5)
%         hold on
        plot(ksT,pdmin_all,'--','LineWidth',1.5)
        % plot(ksT01,pdmin_zErrCut,'--','LineWidth',1.5)
        plot(ksT,pdmin_crop,'--','LineWidth',1.5)
%         plot(ksT,pdmin_zErrCut,'--','LineWidth',1.5) % Remove for simpler plot
%         simplified
        yl = ylim;
        plot([1 1]*atmo{kk}.Tmode,[0 pdMax],'--ok','LineWidth',2)
        plot([1 1]*atmo{kk}.T_halfMax,[0 pdHalfMax],'--o','Color',[0.5 0.5 0.5],'LineWidth',2)
        ylim(yl)
        xlim(prctile(rawT(:),[1 95]))
        xlabel('T_{camera} - T_{profile} [K]')
        ylabel('Probability Density')
        title(sprintf('Total dT probability distributions, with varying filters',atmo{kk}.quickName))
        
        % All the things
%          legend('All Pixels','Filt: t, z, zErr','Row-wise minima, all pixels',...
%             'Row-wise minima, Filt: t, z','Row-wise minima, Filt: t, z, zErr',...
%             sprintf('T_{1/2 Max} = %.2f K',atmo{kk}.T_halfMax),'location','northwest')
        
        % Simplified
        legend('All Pixels','Filt: t, z, zErr','Row-wise minima, all pixels',...
            'Row-wise minima, Filt: t, z',...
            sprintf('Final T_{Mode} = %.2f K',atmo{kk}.Tmode),...
            sprintf('T_{1/2 Max} = %.2f K',atmo{kk}.T_halfMax),'location','northwest')

        % -------- Residuals, Sum of squared error and rms error reporting ----------
        figure(rmsFig)
        
        subplot(1,ncols,kk)
        plot(D.z.*0,D.z,'--k','LineWidth',2)
        hold on
        mo(1) = plot(stats_vs_z_all.Mode-atmo{kk}.Tmode,D.z,'Color',co(1,:),'LineWidth',1.5);
        mo(2) = plot(stats_Tcut.Mode(zFiltFlag)-atmo{kk}.Tmode,Z(zFiltFlag),'Color',co(2,:),'LineWidth',1.8);
        title(sprintf('%s: RMSE_{mode} = %.2f',[atmo{kk}.quickName(1:4) '.' atmo{kk}.quickName(end-3:end)],atmo{kk}.modeTrms))
        xlabel('dT_{mode} [K]')
        ylabel('Height above target [m]')
        legend(mo,{'Unfiltered','Filtered'})
        axis tight
        grid on
        
        rmp(kk)=plot(rmsAx,atmo{kk}.Tssq,Z(~ErrNaN),'LineWidth',2);
        grid(rmsAx,'on')
        hold(rmsAx,'on')
        rml{kk} = sprintf('%s:\tRMSE_{all} = %.2f',[atmo{kk}.quickName(1:4) '.' atmo{kk}.quickName(end-3:end)],atmo{kk}.Trms);
        set(rmsAx,'FontSize',12)
%         plot(errax(2),zErrPct,Tcube.z,zUncPct,Tcube.z,'--','Color',co(kk,:))
%         grid on
       
        % Accumulate the useful stats
        fitStats(kk).inputFrames           = tI;
        fitStats(kk).imgROI                = imgROI;
%         fitStats(kk).t                     = reft;
        fitStats(kk).z                     = D.z;
        fitStats(kk).zFilt                 = Z;       
%         fitStats(kk).Mode_vs_t             = stats_v_t.Mode;
%         fitStats(kk).apparentMin_vs_t      = ym;
%         fitStats(kk).rawT_Mode_vs_z_bestFrames_masked_zErrFiltered   = pdstats_Traw_filtered.Mode;
%         fitStats(kk).medianShift           = median(atmo{kk}.Tinterp-pdstats_Traw_filtered.Mode);
%         fitStats(kk).Mode_vs_z_allFrames   = stats_v_z.Mode;
%         fitStats(kk).apparentMin_vs_z_all  = statsMin_v_z.Mode;
        fitStats(kk).nPix_filtered          = zNpixFilt;
        fitStats(kk).pctPix_filter          = zAllPct;
        fitStats(kk).nPix_above100          = ~zFiltFlag;              
        fitStats(kk).Mode_vs_z_raw          = stats_vs_z_all.Mode;
        fitStats(kk).Mode_vs_z_filtered     = stats_Tcut.Mode;
        fitStats(kk).pdfT                   = ksT;
        fitStats(kk).pdf_allPixels          = pdall;
        fitStats(kk).pdf_filteredPixels     = pd_zErrCut_tcut;
        fitStats(kk).pdf_apparentMinAll     = pdmin_all;
        fitStats(kk).pdf_apparentMinFilt    = pdmin_crop;
%         fitStats(kk).apparentMin_vs_z_best  = statsMin_v_z_tcut.Mode;

    end
            % ------------------ Height Error Summary ---------------------
        figure('name','Height Error estimation','position',[50 50 800 1000])
        subplot(2,1,1)
        pcolor(tt,Z,squeeze(zErrMax(:,2,:)))
        shading flat
        colormap(copper(200))
        hold on
        elvls = 50:50:(round(max(max(squeeze(zErrMax(:,2,:)))))-50);
        [C,h]=contour(tt,Z,squeeze(zErrMax(:,2,:)),elvls,'w');
        clabel(C,h,'Color','w');
        colorbar
        title('Estimated height error at plume centerline due to projection geometry [m]')
        xlabel('time [s]')
        ylabel('height [m above target]')

        subplot(2,1,2)
        pcolor(tt,Z,zErrUncertainty)
        shading flat
        colormap(copper(200))
        hold on
        elvls = 50:50:(round(max(zErrUncertainty(:)))-50);
        [C,h]=contour(tt,Z,zErrUncertainty,elvls,'w');
        clabel(C,h,'Color','w');
        colorbar
        title('Estimated eight error uncertainty due to plume distance [m]')
        xlabel('time [s]')
        ylabel('height [m above target]')        

    figure(rmsFig)
    xlabel(rmsAx,'Sum squared error')
    ylabel(rmsAx,'z [m above target]')
    title(rmsAx,'Filtered data fit for all pixels')
    legend(rmp,rml)
    
    if length(atmo)==1
        D.atmo = atmo{1};
    else
        D.atmo = atmo;
    end
    
    if ~isempty(atmoFile)
        fprintf('Updating saved atmospheric profile:\n\t%s\n',atmoFile)
        save(atmoFile,'atmo','fitStats','-append')
    end
end