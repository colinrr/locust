function varargout=getAtmoProfile(satProfiles,roi,timezoneOffset,Tcube,Tidx,Tthresh,outputDir)
% Usage cases:
% atmo = getAtmoProfile(satProfiles,roi,<timezoneOffset>)
%       --> Get atmospheric profiles from hdf file. If only one hdf file is
%       entered for satProfiles, atmo will be a struct for that profile. If
%       more than one hdf path is entered, output will be a cell of profile
%       structures.
%
% [atmo,D] = getAtmoProfile(satProfiles,roi,<timezoneOffset>,<Tcube>,<Tidx>)
%           --> Get atmospheric profiles, and remove them from thermal 
%           image data to get an initial estimate the best thermal data to
%           use for performing the final fit in "fitAtmoProfile"
%
% [atmo,D] = getAtmoProfile(...<outputDir>)
%           --> Save retrieved profiles to the designated directory
%           (filename is generated automatically)
%
% Profiles from the hdf files will be selected based on spatial proximity 
% to target/ROI and temporal proximity to timestamp. Currently calibrated for
% MODIS "MOD" and "MYD", and AIRS "RetStd" and "RetSup". Distinguishes data
% type based on file name, so don't rename the data files...
%           
% INPUT:
%   satProfiles = {nx1} cell, where entries are the path(s) to hdf files
%
%   roi         = array, [1x2 or 1x4], lat/lon limits of interest.
%                  --> IF [1x2]: [lat lon]. Retrieves the single profile 
%                           from each file that is closest to this point
%                  --> IF [1x4]: [lat1 lat2 lon1 lon2]. Retrieves all
%                  profiles from within these bounds
%
% OPTIONAL:
%   timezoneOffset = for adjusting profile UTC time to local
%
%   Tcube           = path to thermal data cube file, or data cube struct.
%                   This is used to perform an initial atmospheric profile
%                   removal on the thermal data, and is useful for setting
%                   constraints and thresholds in "fitAtmoProfile".
%                   It is recommended to use "getThermCube" to create a 
%                   second data cube containing only reference frames of
%                   interest.
%
%   Tidx        = list of subscripts in third dimension of D.T, identifying
%                 which frames ot use for estimating background profile.
%                   --> Default is every 10 frames out of the last 100,
%                   truncated if there are less than 100 frames in D.T,
%                   though it is not recommended to use this setting.
%
%   outputDir   = directory to save output atmospheric profile. If empty,
%               the profiles struct will not be saved to disk.
%
% OUTPUT:
%   atmo    = struct, or cell of structs, with fields dependent on input
%           files. See SDS fields in the nested function 'retrieveProfile'.
% 
%   quickStats = returns a struct with some simple stats from initial fit
%           estimate (requires that thermal data was input).
%
%   opath   = IF outputDir is specified, the atmo profile struct is 
%           automatically saved in that directory, and the full path is
%           output here.
%
%
% C Rowell, May 2020
fprintf('\n========= Get Atmospheric Profile(s) =========\n')


narginchk(2,7)
if nargin<7
    outputDir = [];
end
if nargin<6
    Tthresh = [];
end
if nargin<5
     Tidx = [];
end
if nargin<4
    Tcube = [];
end
if nargin<3
    timezoneOffset = [];
end

%% Get atmo data from satellite hdf files
atmo = cell(size(satProfiles,1),1);
for ff = 1:size(satProfiles,1)
    
    [atmo{ff}] = retrieveProfile(satProfiles{ff},roi);
end

%% Get thermal background profiles from reference images

if ~isempty(Tcube)
    useThermCube = true;

    
    disp('Loading thermal data cube...') % DO I NEED THIS?
    if ischar(Tcube)
        load(Tcube)
        Tcube = D; clear D;
    elseif ~isstruct(Tcube)
        error('Input ''thermFile'' must be of type ''char'' or ''struct''.')
    end

    if isempty(Tidx)
        Tidx = 1:length(Tcube.t);
    end

    if isvector(Tidx)
        refT = Tcube.T(:,:,Tidx);
        refMask = Tcube.mask(:,:,Tidx);
        reft    = Tcube.t(Tidx);
    end
    

    % Get z error
    disp('Height error estimation...')
    [zErr,zErrMax] = zErrorEstimation(Tcube.x,Tcube.z+(Tcube.geom.Ztarg-Tcube.geom.Z0),refMask,Tcube.geom,'bwdist');
    zErrUncertainty = squeeze(zErrMax(:,3,:)-zErrMax(:,1,:));    

    zASL = Tcube.z+Tcube.z0ASL; % Image height vector ASL


    % Cut out pixels in near proximity to plume boundary 
    % (avoid background/transparent pixels)
    clippedMask = refMask;
    if ~isempty(Tthresh)
        clippedMask(refT<Tthresh) = false;
    end
    for ll=1:size(refT,3)
        clippedMask(:,:,ll) = bwdist(~refMask(:,:,ll)) > 20;
    end
    
    

    fprintf('Processing thermal data for %i atmospheric profile(s)...\n',length(atmo))
    for kk=1:length(atmo)
        atmo{kk}.Tinterp = interp1(atmo{kk}.Height(~isnan(atmo{kk}.Height)),...
            atmo{kk}.Temperature(~isnan(atmo{kk}.Temperature)),zASL,'spline');
        
        clippedT = (refT - atmo{kk}.Tinterp);
        clippedT(clippedMask==0) = NaN; % Subtract out profile and apply mask
        
        % Temperature vector for PDF
        varioCamTres = 0.05; % VarioCam temperature resolution
        ksT = nanmin(clippedT(:)):varioCamTres*4:nanmax(clippedT(:));

        % Minima across all horizontals
        dTmin = squeeze(nanmin(clippedT,[],2));
        % Distributions for all height over time
        [pd_vs_t,~,pdall,stats_v_t] = ksdensityND(clippedT,ksT,[1 2]); % All
        % Probability distrubution of minima over time
        [pdmin_vs_t,~,pdmin,statsMin_v_t] = ksdensityND(dTmin,ksT,1); 
        
        % Distributions for all time over height
        [pd_vs_z,~,~,stats_v_z] = ksdensityND(clippedT,ksT,[2 3]);
        [pdmin_vs_z,~,~,statsMin_v_z] = ksdensityND(dTmin,ksT,2); 
        
%        run sandbox_scripts/atmoProfile_sandboxing
        % ---------- exponential fit vs time? --------------
         %-> ALL T
        x = reft;
        y = double(stats_v_t.Mode-min(stats_v_t.Mode));
        [f0,gof] = fit(x,y,'exp2');
        f0y = f0(x)+double(min(stats_v_t.Mode));
         %-> ALL T
        ym = double(statsMin_v_t.Mode); %-min(statsMin_v_t.Mode));
%         f1 = fit(x,ym,'exp2');         
            
        % ---- Quick estimate of "better" regions for time/height filter ----
        % -> for visualization/data QC only
            GRADprctile = 40;
            SIGprctile  = 75;
             % -> Small gradient/std. dev for mode vs time
            gradMode = abs(gradient(smooth(y,max([25 round(length(y)/5)])),mean(diff(reft))));
            dModePrct = prctile(gradMode,GRADprctile);
            ModeSigPrct = prctile(stats_v_t.StdDev,SIGprctile);
            qct_flag = and(gradMode<dModePrct, stats_v_t.StdDev<ModeSigPrct);
            
            % Cut out some frames with the lowest Mode gradients and standard deviations
            T_tcut = clippedT(:,:,qct_flag);
            [pd_vs_z_tcut,~,~,Stats_tcut] = ksdensityND(T_tcut,ksT,[2 3]);
            [~,~,~,statsMin_v_z_tcut] = ksdensityND(dTmin(:,qct_flag),ksT,2); 
            TstdLims = [nanmean(Stats_tcut.Mode)-5*nanmean(Stats_tcut.StdDev); nanmean(Stats_tcut.Mode)+5*nanmean(Stats_tcut.StdDev)];
            Tcut_in = reft(qct_flag);
            
            % Approximate height difference for 1.5 K uncertainty
            TuncLevel = abs([1.5]/mean(gradient(atmo{kk}.Tinterp,mean(diff(Tcube.z)))));
            T_tcut_zErr = T_tcut;
            T_tcut_zErr(zErr(:,:,qct_flag)>TuncLevel) = NaN;
            [pd_vs_z_tcut_zErr,~,~,Stats_tcut_zErr] = ksdensityND(T_tcut_zErr,ksT,[2 3]);
            
            Traw_filtered = refT(:,:,qct_flag).*refMask(:,:,qct_flag);
            Traw_filtered(Traw_filtered==0) = NaN;
            Traw_filtered(zErr(:,:,qct_flag)>TuncLevel) = NaN;
            [pdraw,~,~,pdstats_Traw_filtered] = ksdensityND(Traw_filtered,nanmin(Traw_filtered(:)):varioCamTres*4:nanmax(Traw_filtered(:)),[2 3]);
            % Get
             % -> Similar for height?
%              dModePrc
             % -> zError stats

    %% PLOTTING PDFS/Z Error stats FOR EACH PROFILE
    
        % Plot total frame dT distrubution over time
        co = get(gca,'ColorOrder');
        figure('name',sprintf('Probability Density: T - %s',atmo{kk}.quickName),'position',[20 50 1000 1000])
        tightSubplot(2,1,1,[],0.07)
        pcolor(reft,ksT,pd_vs_t./max(pd_vs_t,[],1))
        shading flat
        colormap(gray(150))
        title('PDF''s of thermal frames (all height) versus time, atmospheric profile removed')
        hold on
        % Rectangle limits for "best" data
        Rt = plot([reft(find(qct_flag,1,'first')).*[1;1],reft(find(qct_flag,1,'last')).*[1;1]],TstdLims,'c','LineWidth',2);
        plot([reft(qct_flag) reft(qct_flag)],repmat(TstdLims',[sum(qct_flag) 1]),'.-c')
        % Mode values and fit curve
        pdfp(1:2) = plot(x,stats_v_t.Mode,'.',x,f0y,'Color',co(1,:),'LineWidth',1.5);
        % plot(Tcube.t,dTminzI,'r')
        ylabel('dT [K]')
        % Apparent minima values
        pdfp(3) = plot(x,ym,'.','Color',co(2,:),'LineWidth',1.5);
        xlabel('Time [s]')
        legend([Rt(1) pdfp],{'Low-variation data','PDF Mode',sprintf('Mode exp. fit, R^2=%.2f',gof.rsquare),'Apparent minima'})
        set(gca,'FontSize',12)
        
            % Plot total (all time) distrubtion vs height
        tightSubplot(2,3,4,0.1,0.07)
        pcolor(ksT,Tcube.z,(pd_vs_z./max(pd_vs_z,[],1))')
        shading flat
        hold on
        plot(stats_v_z.Mode,Tcube.z,'.','Color',co(1,:))
        plot(statsMin_v_z.Mode,Tcube.z,'.','Color',co(2,:))
        xlabel('dT [K]')
        ylabel('z [m above vent]')
        xlim(prctile(clippedT(:),[2 95]))
        title('Normalized T (all time) PDF vs height')
        set(gca,'FontSize',12)
        
        tightSubplot(2,3,5,0.1,0.07)
        pcolor(ksT,Tcube.z,(pd_vs_z_tcut./max(pd_vs_z_tcut,[],1))')
        shading flat
        hold on
        plot(Stats_tcut.Mode,Tcube.z,'.','Color',co(1,:))
        plot(statsMin_v_z_tcut.Mode,Tcube.z,'.','Color',co(2,:))
        xlabel('dT [K]')
        ylabel('z [m above vent]')
        xlim(prctile(clippedT(:),[2 95]))
        title('Normalized T (low-var) PDF vs height')
        set(gca,'FontSize',12)

        tightSubplot(2,3,6,0.1,0.07)
        pcolor(ksT,Tcube.z,(pd_vs_z_tcut_zErr./max(pd_vs_z_tcut_zErr,[],1))')
        shading flat
        hold on
        plot(Stats_tcut_zErr.Mode,Tcube.z,'.','Color',co(1,:))
%         plot(statsMin_v_z_tcut.Mode,Tcube.z,'.','Color',co(2,:))
        xlabel('dT [K]')
        ylabel('z [m above vent]')
        xlim(prctile(clippedT(:),[2 95]))
        title('Normalized T (low-var, z Error filter) PDF vs height')
        set(gca,'FontSize',12)

        % Record a few stats
        initialStats(kk).inputFrames           = Tidx;
        initialStats(kk).bestFrames            = Tidx(qct_flag);
        initialStats(kk).t                     = reft;
        initialStats(kk).zASL                  = zASL;
        initialStats(kk).Mode_vs_t             = stats_v_t.Mode;
        initialStats(kk).apparentMin_vs_t      = ym;
        initialStats(kk).rawT_Mode_vs_z_bestFrames_masked_zErrFiltered   = pdstats_Traw_filtered.Mode;
        initialStats(kk).medianShift           = median(atmo{kk}.Tinterp-pdstats_Traw_filtered.Mode);
        initialStats(kk).Mode_vs_z_allFrames   = stats_v_z.Mode;
        initialStats(kk).apparentMin_vs_z_all  = statsMin_v_z.Mode;        
        initialStats(kk).Mode_vs_z_bestFrames  = Stats_tcut.Mode;
        initialStats(kk).apparentMin_vs_z_best = statsMin_v_z_tcut.Mode;
    end
    
     %Plot RMS ERROR after fit
%     figure(rmsFig)
%     legend(rmp,rml)

    
else
    useThermCube = false;
    initialStats = [];
end

    %% PLOTTING
    
    % ------- Atmo profile and target locations ---------
    % --> Locations
    % --> Times
    figure
    scatter(roi(2),roi(1),40,'r','filled')
    hold on
    text(roi(2),roi(1),'Target','HorizontalAlignment','left')
    for kk=1:length(atmo)
        atmo{kk}.localTime = atmo{kk}.UTC_Time + timezoneOffset/24;
        
        if useThermCube
            time_offset = atmo{kk}.localTime - Tcube.t0_local;
            scatter(atmo{kk}.Longitude,atmo{kk}.Latitude,40,abs(time_offset),'filled')
        else
            co = get(gca,'ColorOrder');
            scatter(atmo{kk}.Longitude,atmo{kk}.Latitude,40,atmo{kk}.localTime,'filled')
        end
        text(atmo{kk}.Longitude,atmo{kk}.Latitude,...
            sprintf('%s\n%s',atmo{kk}.quickName,datestr(atmo{kk}.localTime)),'HorizontalAlignment','left')

    end
    xlabel('Longitude')
    ylabel('Latitude')
    title('Atmospheric profile and target coordinates')
    
if useThermCube

    % Show summary of zError
    figure('name','Height Error estimation','position',[50 50 800 1000])
    subplot(2,1,1)
    pcolor(reft,Tcube.z,squeeze(zErrMax(:,2,:)))
    shading flat
    colormap(copper(200))
    hold on
    elims = round(max(max(squeeze(zErrMax(:,2,:))))).*[0.1 0.9];
    elvls = linspace(elims(1),elims(2),9);
%     elvls = 50:50:(round(max(max(squeeze(zErrMax(:,2,:)))))-50);
    [C,h]=contour(reft,Tcube.z,squeeze(zErrMax(:,2,:)),elvls,'w');
    clabel(C,h,'Color','w');
    colorbar
    title('Estimated Z error at plume centerline due to plume diameter + projection geometry [m]')
    xlabel('time [s]')
    ylabel('height [m above target]')
    
    subplot(2,1,2)
    pcolor(reft,Tcube.z,zErrUncertainty)
    shading flat
    colormap(copper(200))
    hold on
    elims = round(max(zErrUncertainty(:))).*[0.1 0.9];
    elvls = linspace(elims(1),elims(2),9);
%     elvls = 50:50:(round(max(zErrUncertainty(:)))-50);
    [C,h]=contour(reft,Tcube.z,zErrUncertainty,elvls,'w');
    clabel(C,h,'Color','w');
    colorbar
    title('Estimated Z error uncertainty due to plume-camera distance [m]')
    xlabel('time [s]')
    ylabel('height [m above target]')
    
    % Plot atmo profiles with modal profiles, shifted by Median diff
    figure
    count = 1;
    for kk=1:length(atmo)
        pl(count) = plot(atmo{kk}.Temperature,atmo{kk}.Height,'o--','LineWidth',1,'Color',co(kk,:));
        hold on
        pl(count+1) = plot(atmo{kk}.Tinterp,zASL,'LineWidth',2,'Color',co(kk,:));
        plab{count} = ['Raw ' atmo{kk}.quickName];
        plab{count+1} = ['Interp ' atmo{kk}.quickName];
        pl(count+2) = plot(pdstats_Traw_filtered.Mode+median(atmo{kk}.Tinterp-pdstats_Traw_filtered.Mode),zASL,'Color',co(kk,:)*0.45,'LineWidth',2);
        plab{count+2} = 'Mode, raw T, median shift';
        
        count = count+3;
    end
    pl(count) = plot(pdstats_Traw_filtered.Mode,zASL,'k','LineWidth',2);
    plab{count} = 'Mode, raw T, filtered';
    title('Atmospheric profiles')
    xlabel('T [K]')
    ylabel('z [m ASL]')
    ylim([min(zASL)-1000 max(zASL)+1000])
    legend(pl,plab)    
    
 
else % No input thermal cube
    % Just plot atmo profiles
    figure
    count = 1;
    for kk=1:length(atmo)
        pl(count) = plot(atmo{kk}.Temperature,atmo{kk}.Height,'o--','LineWidth',1,'Color',co(kk,:));
        plab{count} = ['Raw ' atmo{kk}.quickName];
        hold on
        if useThermCube
            pl(count+1) = plot(atmo{kk}.Tinterp,zASL,'LineWidth',2,'Color',co(kk,:));
            plab{count+1} = ['Interp ' atmo{kk}.quickName];
            count = count+2;
        else
            count = count+1;
        end
    end
%     plot(Tmin_stacked,zASL,'k','LineWidth',2);
    title('Atmospheric profiles')
    xlabel('T [K]')
    ylabel('z [m ASL]')
    if useThermCube; ylim([min(zASL)-1000 max(zASL)+1000])
    else; axis tight; end
    grid on
    legend(pl,plab)

end

% Sort output here...
% quickStats;

if ~isempty(outputDir)
    oname = sprintf('atmoProfiles_x%i_%s.mat',length(atmo),datestr(now,'YY-mm-dd'));
    opath = fullfile(outputDir,oname);
    fprintf('Writing atmospheric profile struct to:\n\t%s\n',opath)
    save(opath,'atmo','initialStats')
end
if nargout>=1
    varargout{1} = atmo;
end
if nargout>=2
    if useThermCube
        varargout{2} = initialStats;
    else
        varargout{2} = [];
    end
end
if nargout==3
    if ~isempty(outputDir)
        varargout{3} = opath;
    else
        varargout{3} = [];
    end
end
end


% ==================================================================
% ++++++++++++++++++++++ END MAIN FUNCTION +++++++++++++++++++++++++
% ==================================================================



%% Retrieve Atmospheric Profiles from hdf
function [Ss] = retrieveProfile(hdfile,roi)
% SS = RetrieveRawMODIS(hdfile,sds_fields,roi)
% IN:   hdfile     = full path to MODIS hdf file
%       sds_fields = cell array, with each entry being a two element cell
%               containing {output field name, field name from HDF file}
%               eg. {'Temperature', 'Retrieved_Temperature_Profile'}
%       roi        = OPTIONAL 2 array, [1x2 or 1x4], lat/lon limits
%               of interest.
%                  --> IF [1x2]: [lat lon]. Retrieves the single profile 
%                           from each file that is closest to this point
%                  --> IF [1x4]: [lat1 lat2 lon1 lon2]. Retrieves all
%                  profiles from within these bounds
%                  --> IF []: empty, retrieves ALL data
%   ---> Doesn't seem to work with hdfread right now...
% OUT:  SS = Data structure, where each entry is a structure containing data and
%            attributes for the named field
%
% C Rowell, September 2018
    fprintf('Retrieving data from file:\n\t%s\n',hdfile)
    
    if nargin<2
        roi = [];
    end    
 %%
    
    % MODIS profile fields - Grabs lat/lon/pressure automatically
sdsMOD = {...
    {'UTC_Time',            'Scan_Start_Time'},
    {'Moisture',            'Retrieved_Moisture_Profile'},
    {'Temperature',         'Retrieved_Temperature_Profile'},
    {'WV_Mixing_Ratio',     'Retrieved_WV_Mixing_Ratio_Profile'}, % g/kg
    {'Height',              'Retrieved_Height_Profile'},
    {'Surface_Pressure',    'Surface_Pressure'},
    {'Surface_Elevation',   'Surface_Elevation'},
    {'Tropopause_Height',   'Tropopause_Height'},
    {'QC',                  'Quality_Assurance'}
    };

% AIRS Standard Profile fields
sdsAStd = {
%     {'Pressure',            'pressStd'},
    {'UTC_Time',            'Time'},
    {'Temperature',         'TAirStd'},
    {'TempErr',             'TAirStdErr'},
    {'Temp_QC',             'TAirStd_QC'},
    {'Height',              'GP_Height'},
    {'Height_QC',              'GP_Height_QC'},
    {'nSurface',            'nSurfStd'},
    {'Surface_Pressure',    'PSurfStd'},
    {'SurfAirTemp',         'TSurfAir'},
    {'SurfAirTempErr',      'TSurfAirErr'},
    {'RelHumidity',         'RelHum'},
    {'RelHumidity_QC',      'RelHum_QC'},
%     {'WV_Mixing_Ratio',     'H2OMMRLevStd'},
%     {'WV_Mixing_Ratio_Err', 'H2OMMRLevStdErr'},
%     {'WV_Mixing_Ratio_QC',  'H2OMMRLevStd_QC'},
    % See H2OMMRSatLevStd for saturation mixing ratios
    {'nBest_QC',            'nBestStd'},
    {'nGood_QC',            'nGoodStd'}
    };

% AIR Supplemental Profile fields
sdsASup = {
%     {'Pressure',            'pressSup'},
    {'UTC_Time',            'Time'},
    {'Temperature',         'TAirSup'},
    {'TempErr',             'TAirSupErr'},
    {'Height',              'GP_HeightSup'},
    {'Height_QC',           'GP_HeightSup_QC'},
    {'nSurface',            'nSurfStd'},
    {'Surface_Pressure',    'PSurfStd'},
    {'SurfAirTemp',         'TSurfAir'},
    {'SurfAirTempErr',      'TSurfAirErr'},
    {'RelHumidity',         'RelHum'},
    {'RelHumidity_QC',      'RelHum_QC'},
    %     {'WV_Mixing_Ratio',     'H2OMMRLevSup'},
%     {'WV_Mixing_Ratio_Err', 'H2OMMRLevSupErr'},
%     {'WV_Mixing_Ratio_QC',  'H2OMMRLevSup_QC'},
    {'nBest_QC',            'nBestSup'},
    {'nGood_QC',            'nGoodSup'}
    };

    [~,fname,~] = fileparts(hdfile);
    switch true
        case or(contains(fname,'MOD'),contains(fname,'MYD'))
            pressure_field  = 'Pressure_Level';
            sds_fields      = sdsMOD;
            correctFields = true;
            quickname = fname(1:22);
%             [profiles{ff}] = retrieveProfile(satProfiles{ff},'Pressure_Level',sdsMOD,roi); %fliplr(atm_roi));

%             [ventN,ventE,~] = ell2utmWGS84(vent(2), vent(1));
%             [atmo,lref] = getClosestProfile(modis,[ventE ventN],'utm');
            
        case contains(fname,'RetStd')
            pressure_field  = 'pressStd';
            sds_fields      = sdsAStd;
            correctFields = false;
            quickname = fname(1:19);
%             [profiles{ff}] = retrieveProfile(satProfiles{ff},'pressStd',sdsAStd,roi);
            
        case contains(fname,'RetSup')
            pressure_field  = 'pressSupp';
            sds_fields      = sdsASup;
            correctFields = false;
            quickname = fname(1:19);
%             [profiles{ff}] = retrieveProfile(satProfiles{ff},'pressSup',sdsASup,roi);
    end
    %%
  
    Lat  = double(hdfread(hdfile,'Latitude'));
    Lon  = double(hdfread(hdfile,'Longitude'));
    P    = hdfread(hdfile, pressure_field);
    LatSizes = size(Lat);
    
    II = hdfinfo(hdfile);
    SDSinfo = II.Vgroup.Vgroup(2).SDS;
    SDSnames = {SDSinfo.Name}';
    
    Vdatinfo = II.Vgroup.Vgroup(2).Vdata;
    Vdatnames = {Vdatinfo.Name}';
    
    GLi = strcmp({II.Vgroup.Vgroup.Name},'Geolocation Fields');
    GLinfo = II.Vgroup.Vgroup(GLi).SDS;
    GLnames = {II.Vgroup.Vgroup(GLi).SDS.Name};
    
    Ss.dataFile = hdfile;
    Ss.quickName = quickname;

    % Retreiving closest point, region, or all?
    switch numel(roi)
        case 0
            getAllProfiles = true;
            Ss.Latitude  = Lat;
            Ss.Longitude = Lon;
        case 2
            getAllProfiles = false;
            [~,~,latlonI] = closest2d(roi(2),roi(1),Lon(:),Lat(:));
            Ss.Latitude  = Lat(latlonI);
            Ss.Longitude = Lon(latlonI);
        case 4
            getAllProfiles = false;
            latCheck = and(Lat>=roi(1),Lat<roi(2));
            lonCheck = and(Lon>=roi(3),Lon<roi(4));
            latlonI = find(and(latCheck,lonCheck));
            if isempty(latlonI)
                error(sprintf('No profiles found in ROI:\n  File:\t%s\n  ROI:\t%s\n',hdfile,mat2str(roi)))
            end
            Ss.Latitude  = Lat(latlonI);
            Ss.Longitude = Lon(latlonI);
    end
    
    if iscell(P)
        P = P{1};
    end
    Ss.Pressure = double(P);
    
    for ss = 1:length(sds_fields)
        % Vdata, Geolocation, or regular SDS?
%         fprintf('%s\n',sds_fields{ss}{2})
        vI = find(strcmp(Vdatnames,sds_fields{ss}{2}));
        gI = find(strcmp(GLnames,sds_fields{ss}{2}));        
        hI = find(strcmp(SDSnames,sds_fields{ss}{2}));

        % Get field meta data
        switch false
            case isempty(hI)
                fieldInfo = SDSinfo(hI);
            case isempty(vI)
                fieldInfo = Vdatinfo(vI);
            case isempty(gI)
                fieldInfo = GLinfo(gI);
        end
        
        % Get data field
        if getAllProfiles
%             for ss = 1:size(sds_fields,1)
%                 fprintf('%s\n',sds_fields{ss}{2})
                Ss.(sds_fields{ss}{1}) = hdfread(hdfile,sds_fields{ss}{2});
%             end

        else

%             for ss = 1:length(sds_fields)
%                 if and(~isempty(gI),isempty(hI))
%                     hI = gI;
%                 end

            [subI,subJ] = ind2sub(size(Lat),latlonI);
            
            % Get indices of the selected data from ROI
            if or(~isempty(hI),~isempty(gI)) % For SDS and Geolocation fields
                nDims = length(fieldInfo.Dims);
                switch nDims
                    case 2
                        hdfIndex = {[min(subI) min(subJ)],[],[numel(subI) numel(subJ)]};
                    case 3
                        dimSizes = [fieldInfo.Dims.Size];
                        vDim = find(~any(dimSizes==LatSizes',1));
                        if numel(vDim)~=1
                            error('Cannot recognize vertical dimnsion.')
                        end

                        switch vDim
                            case 1
                                hdfIndex = {[1 min(subI) min(subJ)],[],[fieldInfo.Dims(vDim).Size numel(subI) numel(subJ)]};
                            case 2
                                hdfIndex = {[min(subI) 1 min(subJ)],[],[numel(subI) fieldInfo.Dims(vDim).Size numel(subJ)]};
                            case 3
                                hdfIndex = {[min(subI) min(subJ) 1],[],[numel(subI) numel(subJ) fieldInfo.Dims(vDim).Size]};
                        end         
                end
    
                % Retrieve data
                Ss.(sds_fields{ss}{1}) = hdfread(hdfile,sds_fields{ss}{2},'Index',hdfIndex);



                elseif ~isempty(vI) % For VData fields
                    Ss.(sds_fields{ss}{1}) = hdfread(hdfile,sds_fields{ss}{2});

                else % Something missing...
                    Ss.(sds_fields{ss}{1}) = [];
                end
%             end
        end
        
        %--- Apply any necessary corrections ----
        
        % Get correction attributes
        if correctFields
            corr_fields = {'valid_range','_FillValue','scale_factor','add_offset'};
        else
            corr_fields = {'_FillValue'};
        end
        
        for kk=1:length(corr_fields)
            [~,li] = ismember(corr_fields{kk},{fieldInfo.Attributes.Name}');
            corrs.(genvarname(corr_fields{kk})) = fieldInfo.Attributes(li).Value;
        end
        
        % Convert to double and clear out fill values
        Ss.(sds_fields{ss}{1}) = double(Ss.(sds_fields{ss}{1}));
%         corrs.valid_range = double(corrs.valid_range);
        nn = Ss.(sds_fields{ss}{1})==corrs.x_FillValue;
        Ss.(sds_fields{ss}{1})(nn) = NaN;
%         ss.x_FillValue = NaN;
        
        if correctFields
            % valid ranges
            % Check if either neither fill value nor inside valide range?
            valid_check = ~or(Ss.(sds_fields{ss}{1})(:)==corrs.x_FillValue,...
                and(Ss.(sds_fields{ss}{1})(:)>corrs.valid_range(1),Ss.(sds_fields{ss}{1})(:)<corrs.valid_range(2)));
            if any(valid_check)
%                 warning(sprintf('Invalid values found in:\t%s\n',sds_fields{ss}{2}))
                Ss.invalid = valid_check;
            end       

            %-> Offsets and scales
            Ss.(sds_fields{ss}{1}) = corrs.scale_factor*(double(Ss.(sds_fields{ss}{1})) - corrs.add_offset);
%             ss.valid_range = ss.scale_factor*(double(ss.valid_range) - ss.add_offset);
        end
        
    end
    
    % Time check
    if ismember('UTC_Time',fieldnames(Ss))
        t = Ss.UTC_Time;
        Ss.UTC_Time = datenum([1993 1 1 0 0 0])+t/86400;
    end
    fprintf('\n')

end

%% closest2d
function [subI,subJ,idx] = closest2d(xi,yi,X,Y)
% [subx,suby,idx] = closest2d(xi,yi,X,Y);
% For gridded or non-gridded 2d data sets, find the point closest to the 
% input value using a pythogorean distance...
% Crude for geographic reference systems, but can use geodetic2are for
% lat/lon coordinates.
% IN:   xi   = reference position(s) x
%       yi   = reference position(s) y
%       X    = vector or grid of x coords
%       Y    = vector or grid of y coords
%
% NOT IMPLEMENTED:
%       crs  = add 'geo' to use the geodetic2aer function for lat/lon coords
%
% OUT:  subI = Dim1 SUBSCRIPT of closest value in X,Y
%       subJ = Dim2 SUBSCRIPT of closest value in X,Y
%       idx  = INDEX of closest value in X and J
%
% C Rowell, July 2018

% if nargin<5
%     crs = '';
% end

% if isempty(crs)

idx = zeros(size(xi));
for ii = 1:size(idx,2)
    for jj  = 1:size(idx,1)
        D = sqrt((xi(ii,jj) - X).^2 + (yi(ii,jj) - Y).^2);
        [~,idx(ii,jj)] = min(D(:));
    end
end

% else
%     [az,elev,D] = geodetic2aer(X,Y,h?,xi,yi,hi?);
% end

[~,idx] = min(D(:));

[subI, subJ] = ind2sub(size(D),idx);

end



%% write atmo
function writeAtmo(S,opath)

ff = fieldnames(S);
heads = {};
name_row = {};
unit_row = {};
dat = [];


for ii = 1:length(ff)
    fn = ff{ii};
    if ~strcmp(fn,'units')
        if ischar(S.(fn))
            oname = sprintf('%s [%s]',fn,S.units.(fn));
            hline = sprintf('%s:\t%s\n',oname,S.(fn));
            heads = [heads; {hline}];
        elseif isscalar(S.(fn))
            oname = sprintf('%s [%s]',fn,S.units.(fn));
            hline = sprintf('%s:\t%f\n',oname,S.(fn));
            heads = [heads; {hline}];
        elseif isvector(S.(fn))
            name_row = [name_row; {fn}];
            unit_row = [unit_row; {S.units.(fn)}];
            dat = [dat S.(fn)];
        end
    end
end

delim = '\t';
fileID = fopen(opath,'w');
for ii = 1:length(heads)
    fprintf(fileID,heads{ii});
end
N = numel(name_row);
nameformat = [repmat('%s\t',[1,N]) '\n'];
fprintf(fileID,nameformat,name_row{:});
fprintf(fileID,nameformat,unit_row{:});

for ii = 1:size(dat,1)
    s = string(dat(ii,:));
    s(ismissing(s)) = 'NaN';
    row_format = [repmat('%f\t',[1,N]) '\n'];
    fprintf(fileID,row_format,s);
end
fclose(fileID);
end


