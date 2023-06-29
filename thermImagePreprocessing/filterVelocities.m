function [V,Vspec] = filterVelocities(V,mask,band,filtlims,npoles,testMode,plotFlag)
%   [Vo,Vspec]  = filterVelocities(V,band,filtlims,npoles,testflag)
%   V           = velocity struct output from thermOpticFlow2O. If entered
%                 as a path to a data file, the output will be written to a
%                 new file in the same directory, IF testMode = FALSE.
%
%   mask        = 3D mask array from thermal cube corresponding to V
%   
%   band        = 'low', 'high','bandpass'; DEFAULT = 'low'
%
%   filtlims    = [Tlo Thi] or [T]: period(s) in seconds for filter limits.
%                 Default is lowpass, with cutoff period = mode([u^2 +
%                 v^2]^1/2)/dx, where dx is grid spacing.
%
%   npoles      = number of poles for butterworth filter
%
%   testMode    = true/false. When True mode, will calculate and plot
%                 spectra, running a test filter on a small subset of data
%                 to see the results. Default TRUE.
%                 
%   plotFlag    = true/false to plot diagnostic spectra and timeseries
%
% If second output specified (Vspec), will output raw and filtered spectra in
% second [2x1] struct

%% SETUP/PARSE
    disp('------------ Filtering Velocity Cube -----------')

    narginchk(2,7)
    if nargin<5
        plotFlag = true;
    end
    if isempty(plotFlag)
        plotFlag = true;
    end
    if nargin<6
        testMode = true;
    end
    if isempty(testMode)
        testMode = true;
    end
    if nargin<5
        npoles = 4;
    end
    if isempty(npoles)
        npoles = 4;
    end
    if nargin<4
        filtlims = [];
    end
    if nargin<3
        band = 'low';
    end
    if isempty(band)
        band = 'low';
    end

    if ischar(V)
        if exist(V,'file')
            fprintf('Loading:\n\t%s\n',V)

            if ~testMode
                fprintf('\tWrite mode ON\n')
                writeMode = true;
                ifile = V;
            else
                writeMode = false;
            end
            load(V)
        end
    elseif ~isstruct(V)
        error('"V" input not recognized')
    else
        writeMode = false;
    end

%% Get filtlims defaults 
    dt = mean(diff(V.t));
    uGrid = V.dz/dt;

        % V/uGrid PDF's and mode
        disp('Calculating velocity probability distribution...')
        Vall = sqrt(V.Vz(mask).^2 + V.Vx(mask).^2);
        uVec = 0:0.1:max([uGrid prctile(Vall,99)]);
        V_pdf = zeros(length(uVec),10);
        for ii=1:size(V_pdf,2)
            Vsamp = randsample(Vall,1e4);
            V_pdf(:,ii) = ksdensity(Vsamp,uVec);
        end
        V_pdf = mean(V_pdf,2);
        [~,vI] = max(V_pdf);
        
        % Modal Courant ratio
        CMode = uVec(vI)/uGrid;
        % Modal pixel traverse time
        pixelTime = 1/CMode*dt;
        
        if isempty(filtlims)
            filtlims = pixelTime;
            band = 'low';
        end

%% Build filter
    Fs = 1./dt;
    nyquist = Fs/2;
    [b,a] = butter(npoles,1./filtlims./nyquist,band);
   
    cubeSz = size(V.Vx);

%% Visualization: Velocity PDF, Sample Spectra + Time series before/after

        % Sample velocity time series for plotting
        disp('Retrieving sample time series...')
        Nspec = 100;
%         allMask = all(mask,3);
        allMask = sum(mask,3)>=0.75*size(mask,3);
        
        sampIdx = randsample(find(allMask),Nspec);
        [sampI,sampJ] = ind2sub(size(allMask),sampIdx);
        u = zeros(size(V.Vz,3),Nspec);
        w = u;
        masksamp = logical(u);
        for ii=1:Nspec
            u(:,ii) = squeeze(V.Vx(sampI(ii),sampJ(ii),:))';
            w(:,ii) = squeeze(V.Vz(sampI(ii),sampJ(ii),:))';
            masksamp(:,ii) = squeeze(mask(sampI(ii),sampJ(ii),:))';
        end
        U = sqrt(u.^2 + w.^2);
        
        % Get pre-filtered Spectra
        disp('Calculating sample spectra...')
%         pxx_raw = zeros(size(u,1),3);
        [pxxu,~] = pmtm(u,[],[],1/dt);
        [pxxw,~] = pmtm(w,[],[],1/dt);
        [pxxU,F] = pmtm(U,[],[],1/dt);
        pxxRaw = [mean(pxxu,2) mean(pxxw,2) mean(pxxU,2)];

    if plotFlag

        % Apply filter to sample time series, get spectra
        uf = filtfilt(b,a,double(u));
        wf = filtfilt(b,a,double(w));
        Uf = filtfilt(b,a,double(U));
        [pxxuf,~] = pmtm(uf,[],[],1/dt);
        [pxxwf,~] = pmtm(wf,[],[],1/dt);
        [pxxUf,~] = pmtm(Uf,[],[],1/dt);
        pxxFilt = [mean(pxxuf,2) mean(pxxwf,2) mean(pxxUf,2)];
        
        % --- PLOTS -----
        % "courant number" distribution
        figure
        plot(uVec/uGrid,V_pdf,'LineWidth',2)
        grid on
        hold on
        yl = ylim;
        cm=plot(CMode*[1 1],yl,'--','Color',[0.3 0.3 0.3],'LineWidth',2);
        xlabel('\textbf{$|\vec{u}|\frac{dt}{dz}$}','Interpreter','latex')
        ylabel('PDF')
        title('Ratio: Optic Flow Velocities to image grid speed')
        legend(cm,{sprintf('Ratio mode = %.2f',CMode)})
        set(gca,'FontSize',16)
        
        % spectra
        figure('position',[200 100 1200 900])
        co = get(gca,'ColorOrder');

        loglog(F,pxxRaw,'LineWidth',2)
        set(gca,'ColorOrder',[co(1:size(pxxRaw,2),:)*0.6;co(1:size(pxxRaw,2),:)])
        hold on
%         loglog(F,wAll,'LineWidth',1.5)
%         loglog(F,Uall,'LineWidth',1.5)
        axis tight;
        yl = ylim;
        
%         set(gca,'ColorOrder',co,'ColorOrderIndex',1)
        loglog(F,pxxFilt,'--','LineWidth',1.6)
        plot(1/pixelTime*[1 1],yl,'-.','Color',[0.3 0.3 0.3],'LineWidth',2)
        ylim(yl)
        grid on
        xlabel('Frequency [Hz]')
        ylabel('Velocity Power Spectral Density')
        if any(pixelTime~=filtlims)
            plot(1/filtlims*[1 1],yl,':k','LineWidth',2)
            legend({'$u_x$','$u_z$','$|\vec{u}|$',...
                    '$u_x$: filtered','$u_z$: filtered','$|\vec{u}|$: filtered',...
                ['$\frac{|\vec{u}_{mode}|}{u_{grid}}dt$ = ' sprintf('%.2f Hz',round(1/pixelTime,2))],sprintf('Cutoff = %s Hz',mat2str(round(1./filtlims,2)))},'Interpreter','latex')
        else
            legend({'$u_x$','$u_z$','$|\vec{u}|$',...
                    '$u_x$: filtered','$u_z$: filtered','$|\vec{u}|$: filtered',...
                ['$\frac{|\vec{u}_{mode}|}{u_{grid}}dt$ = ' sprintf('%.2f Hz',round(1/pixelTime,2))]},'Interpreter','latex')
        end
        set(gca,'FontSize',16)
        
        % Random time series for QC
        Nts = 6;
        tsI = randsample(Nspec,Nts);
        figure('position',[100 100 1000 1000])
        for sp=1:Nts
            tightSubplot(Nts,1,sp)
            plot(V.t,[w(:,tsI(sp)) wf(:,tsI(sp))])
            hold on
            plot(V.t,wf(:,tsI(sp)).*masksamp(:,tsI(sp)),'--','Color',co(2,:))
            axis tight
            grid on
            if sp==1
                legend('Raw','Filtered','Filt + mask')
                title('w: Randomly sampled time series')
            elseif sp~=Nts
                set(gca,'XTickLabel',[])
            elseif sp==Nts
                xlabel('Time [s]')
            end
            ylabel('u_z [m/s]')
        end
        

    end
        

    %% Apply filter to full cube
        Vspec.band        = band;
        Vspec.order       = npoles;
        Vspec.courant     = CMode;
        Vspec.defCutoffHz = 1./pixelTime;       
        Vspec.filtlimsHz  = 1./filtlims;
        Vspec.pxxRaw_uvU  = pxxRaw;
        Vspec.pxxFilt_uvU = pxxFilt;
        Vspec.FreqHz      = F;

    if ~testMode
        disp('Filtering all velocities...')
        if plotFlag
            rawRise = squeeze(max(abs(V.Vz).*mask,[],2));
        end
        V.Vx = reshape(filtfilt(b,a,reshape(double(V.Vx),[prod(cubeSz(1:2)) cubeSz(3)])')',cubeSz);
        V.Vz = reshape(filtfilt(b,a,reshape(double(V.Vz),[prod(cubeSz(1:2)) cubeSz(3)])')',cubeSz);
        
        
        V.filter = Vspec;
        
        if plotFlag
            filtRise = squeeze(max(abs(V.Vz).*mask,[],2));
            
            figure('position',[150 150 1000 1000])
            subplot(2,1,1)
            imagesc(V.t,V.z,rawRise)
            colormap(CubeHelix(200))
            set(gca,'YDir','normal');
            caxis(prctile(abs(rawRise(:)),99)*[0 1])
            set(gca,'XTickLabel',[])
            ylabel('z [m]')
            cb = colorbar;
            cb.Label.String = 'u_z [m/s]';
            title('Velocity Rise Diagram, raw w_{max}')

            tightSubplot(2,1,2)
            imagesc(V.t,V.z,filtRise)
    %         colormap(redblue(200))
            set(gca,'YDir','normal');
            caxis(prctile(abs(rawRise(:)),99)*[0 1])
            xlabel('Time [s]')
            ylabel('z [m]')
            cb = colorbar;
            cb.Label.String = 'u_z [m/s]';
            title('Velocity Rise Diagram, filtered w_{max}')
        end
        

    end

    %% WRITE OUTPUT
    if writeMode
        [odir,iname,ext] = fileparts(ifile);
        oname = [iname '_filtered' ext];
        ofile = fullfile(odir,oname);
        fprintf('Writing filtered velocity struct:\n\t%s\n',ofile)
        save(ofile,'V','Vspec','-v7.3')
    end


end
