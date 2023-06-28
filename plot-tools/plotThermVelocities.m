function varargout = plotThermVelocities(x,z,Vx,Vz,varargin) %,Vmax,opts.idx,pt,opts.roi,opts.plotmode,qcolor)
% plotThermVelocities(x,z,Vx,Vz,varargin)
% OPTIONAL NAME/VALUE PAIR (or Struct) INPUT:
% axis      = specify an axis
% mask      = mask for object of interest- Plots the outline. Leave empty to
%       skip. Either same size as Vx or [size(Vx,1) size(Vx,2) length(opts.idx)]
% trajectory = tracked pixel trajectories or other lines in X,Y,t space.
%               '-> Struct output from "computeTrajectories"
% thermal   = Input thermal data cube - plots thermal values instead of
%               velocity colours, and set velocity to 'vector' mode. Must be
%               identical size to Vx,Vz
% tracks    = struct output from pulseTracker to outline tracked clusters.
% Vmax      = max velocity value in colorbar
% idx       = list of indices in 3rd dimension of velocity cubes
% fdt       = [0.1] pause time between frames (s)
% ROI       = zoom axes limits [y1 y2 x1 x2]
% plotmode  = ['both'] 'colorplot', 'vector', 'both'
% qcolor    = [black, or cyan in thermal mode] 
%             color of quiver vectors ('both' or 'vector' modes)
% Trange    = color axis range for plotting thermal (only relevant w/
%               'Thermal' input)
% qspace = [10] pixel spacing for quiver field
% Qscale    = [2] quiver scale factor
% time      = time vector corresponding to the frames for labeling
%             [length(idx) x 1] or [size(Vx,3) x 1]
% disableV  = true/false; True only valid if 'thermal' option is used.
%             Plots only thermal data with no velocity vectors.
% atmo      = atmospheric profile struct from fitAtmoProfile. Will be
%             removed from any thermal plot, with threshold applied and 
%             mask adjusted.
% colormap  = custom colormap for thermal plot
% fs        = font size (default = 11)
% lw        = line width (default = 1.5)
%
% OPTIONAL OUTPUT: {1}: [1 x 2] axes handles:  [main_axis legend_axis]
%                  {2}: Handle for thermal surface/image object
%                  {3}: Handles for plotted tracks
%%
nargoutchk(0,3)

defQcolor   = [];
defPlotMode = 'both';
defROI      = [];
defPt       = 0.1;
defIdx      = []; 
defVmax     = [];
defTraj     = [];
defmask     = [];
defTherm    = [];
defDpx      = 10; % Quiver decimation factor
defQsc      = 2;  % Quiver scale factor
defTrange   = [240 400];
deft        = [];
defColormap = plasmagrey(150);

defFS = 11; % Default font size
defLW = 1.5; % Default line size

p = inputParser;
addParameter(p,'axis',gca)
addParameter(p,'mask',defmask)
addParameter(p,'trajectory',defTraj)
addParameter(p,'thermal',defTherm)
addParameter(p,'tracks',[])
addParameter(p,'vmax',defmask)
addParameter(p,'idx',defIdx)
addParameter(p,'fdt',defPt)
addParameter(p,'roi',defROI)
addParameter(p,'plotmode',defPlotMode)
addParameter(p,'qcolor',defQcolor)
addParameter(p,'qscale',defQsc)
addParameter(p,'qspace',defDpx)
addParameter(p,'Trange',defTrange)
addParameter(p,'time',deft)
addParameter(p,'disableV',false)
addParameter(p,'colormap',defColormap)
addParameter(p,'colorOrder',[])
addParameter(p,'atmo',[])
addParameter(p,'fs',defFS)
addParameter(p,'lw',defLW)
parse(p,varargin{:})

opts = p.Results;

if isempty(opts.idx)
    opts.idx = 1:size(Vx,3);
end

% Disable if no thermal is input
if opts.disableV && isempty(opts.thermal)
    opts.disableV = false;
elseif and(isempty(Vx),isempty(Vz)) && ~isempty(opts.thermal) && ~opts.disableV
    opts.disableV = true;
end
if ~opts.disableV
    iSz = size(Vx,1);
    jSz = size(Vx,2);
    kSz = size(Vx,3);
else
    iSz = size(opts.thermal,1);
    jSz = size(opts.thermal,2);
    kSz = size(opts.thermal,3);
end

% Thermal input prep
if ~isempty(opts.thermal)
    opts.plotmode = 'vector'; % Overwrites other modes
    ThermLengthCheckI = size(opts.thermal,3)==length(opts.idx);
    ThermLengthCheckV = size(opts.thermal,3)==kSz; 
    assert(and(all([size(opts.thermal,1)==iSz size(opts.thermal,2)==jSz]), ...
        or(ThermLengthCheckI,ThermLengthCheckV)),'Check input thermal dimensions.')    
    if ThermLengthCheckV
        opts.thermal = opts.thermal(:,:,opts.idx);
    end
    
    % Removing atmospheric profile
    if ~isempty(opts.atmo)
        opts.thermal = opts.thermal - opts.atmo.Tinterp - opts.atmo.T_halfMax;
    end
end

% Mask outline setup
if ~isempty(opts.mask)
    % Check dimensions
    MaskLengthCheckI = size(opts.mask,3)==length(opts.idx);
    MaskLengthCheckV = size(opts.mask,3)==kSz;
    assert(and(all([size(opts.mask,1)==iSz size(opts.mask,2)==jSz]), ...
        or(MaskLengthCheckI,MaskLengthCheckV)),'Check input mask dimensions')
    % Get correct frames
    if MaskLengthCheckV
        opts.mask = opts.mask(:,:,opts.idx);
    end
    
    if ~isempty(opts.atmo)
        threshMask = opts.thermal>= -2*(opts.atmo.Tmode-opts.atmo.T_halfMax);
        opts.mask = and(opts.mask,threshMask);
    end
    
    % Get polygons
    for mi = size(opts.mask,3):-1:1
        mpol = mask2poly(opts.mask(:,:,mi));
        if length(mpol)>1
            [~,pi] = max([mpol.Length]);
            mpoly(mi) = mpol(pi);
        elseif ~isempty(mpol)
            try
                mpoly(mi) = mpol;
            catch
                mpol,mi
                pause(1)
            end
        elseif isempty(mpol)
            mpoly(mi) = struct('Length',[],'X',[],'Y',[],'IsFilled',[]);
        end
    end
else
    mpoly = [];
end


if isempty(opts.plotmode)
    opts.plotmode = 'both';
end



assert(all(size(Vx)==size(Vz)),'Vx and Vz must be the same size')

if ~isempty(opts.time)
    if length(opts.time)==kSz
        opts.time = opts.time(opts.idx);
    elseif length(opts.time)~=length(opts.idx)
        error('Length of time vector must correspond to either 3rd dimension of input velocity data or the length of idx.')
    end
end

% Crop velocity cube
if ~opts.disableV
    Vx = Vx(:,:,opts.idx);
    Vz = Vz(:,:,opts.idx);
    if isempty(opts.vmax)
        opts.vmax = max( sqrt(Vx(:).^2 + Vz(:).^2) );
    end
end

if isempty(x)
    x = 1:jSz;
end
if isempty(z)
    z = 1:iSz;
end
if isempty(opts.qcolor)
    if isempty(opts.thermal)
        opts.qcolor = [0.2 0.2 0.2];
    else
        opts.qcolor = [0 0.6 0.6];
    end
end

N = numel(opts.idx);

%% Make a legend
if or( strcmp(opts.plotmode,'colorplot') , strcmp(opts.plotmode,'both') )
    xl = repmat(linspace(-1,1,81),[81 1]);
    yl = xl';
    leg = computeColor(xl,yl);
end

%% Plots

for ii=1:N
    if ii==1 
        ax(1) = opts.axis;
%         axes(ax(1))
    else
%         axes(ax(1))
        cla(ax(1))
    end

    % Plot colorfield and color legend
    if or( strcmp(opts.plotmode,'colorplot') , strcmp(opts.plotmode,'both') )
        img = computeColor(Vx(:,:,ii)./opts.vmax,Vz(:,:,ii)./opts.vmax);
        imagesc(ax(1),x,z,img)
        caxis([-1 1])
        colormap(ax(1),redblue(150))
        set(gca,'YDir','normal')
        set(ax(1),'FontSize',opts.fs)
        xlabel(ax(1),'X [m]')
        ylabel(ax(1),'Z [m]')
        title(ax(1),sprintf('Index = %i',opts.idx(ii)))
        hold(ax(1),'on')
        if ii==1
            axpos  = get(ax(1),'position');
            ax(2)=axes('position',[ axpos(1)+axpos(3)*0.8 axpos(2)+axpos(4)*0.8 axpos(3)*0.2  axpos(4)*0.2 ]);
%         else
%             axes(ax(2))
        end
        imagesc(ax(2),xl(1,:)*opts.vmax,yl(:,1)*opts.vmax,leg)
        set(ax(2),'YDir','normal')
        xlabel('V_x')
        ylabel('V_z')
        if opts.vmax<1
            Vmaxax = fix(Vmax*10)/10;
        else
            Vmaxax = fix(opts.vmax);
        end
        set(ax(2),'XTick',[-Vmaxax Vmaxax])
        set(ax(2),'YTick',[-Vmaxax Vmaxax])
        set(ax(2),'FontSize',opts.fs,'FontWeight','Bold')
        hold(ax(2),'on')
        plot(round([-opts.vmax opts.vmax]),[0 0],'k')
        plot([0 0],round([-opts.vmax opts.vmax]),'k')
%         axes(ax(1))
        daspect(ax(2),[1 1 1])
    end
    
    % Plot thermal?
    if ~isempty(opts.thermal)
%         ss=pcolor(ax(1),x,z,opts.thermal(:,:,ii));
        ss=surf(ax(1),x,z,opts.thermal(:,:,ii),'EdgeAlpha',0);
        TM = max(max(opts.thermal(:,:,ii)));
        
        view([0 0 1])
%         shading(ax(1),'flat')
        colormap(ax(1),opts.colormap)
        xlabel(ax(1),'X [m]')
        ylabel(ax(1),'Z [m]')
        cb=colorbar(ax(1),'location','north');
        cb.Color = [0.9 0.9 0.9];
        if isempty(opts.atmo)
            Tlabel = 'Brightness Temp. [K]';
        else
            Tlabel = 'Excess Temp. [K]';
        end
        cb.Label.String = Tlabel;
        cb.FontSize = opts.fs;
        set(ax(1),'FontSize',opts.fs)
        caxis(ax(1),opts.Trange)
        hold(ax(1),'on')
    else
        TM = 1;
    end

    if isempty(opts.time)
        title(ax(1),sprintf('Index = %i',opts.idx(ii)))
    else
        title(ax(1),sprintf('Index = %i, t = %.2f s',opts.idx(ii),opts.time(ii)))
    end

    daspect(ax(1),[1 1 1])
    
%     if strcmp(opts.plotmode,'both'); hold on; end
    
    % Plot quiver?
    if and( or( strcmp(opts.plotmode,'vector') , strcmp(opts.plotmode,'both') ), ~opts.disableV)
        zv = opts.qspace:opts.qspace:numel(z);
        xv = opts.qspace:opts.qspace:numel(x);
%         quiver(ax(1),x(xv),z(zv),Vx(zv,xv,ii),Vz(zv,xv,ii),opts.qscale,'LineWidth',opts.lw*0.65,'Color',opts.qcolor)
        quiver3(ax(1),x(xv),z(zv),TM.*ones(size(Vx(zv,xv,ii))),Vx(zv,xv,ii),Vz(zv,xv,ii),zeros(size(Vz(zv,xv,ii))),opts.qscale,'LineWidth',opts.lw*0.65,'Color',opts.qcolor)
        hold(ax(1),'on')
    end

        % Plot mask outline?
    if ~isempty(mpoly)
%         hold on
        plot3(ax(1),x(mpoly(ii).X),z(mpoly(ii).Y),TM+mpoly(ii).Y*0,'y','LineWidth',opts.lw)
    end
    
   % Plot trajectories?
    if ~isempty(opts.trajectory)
        tix = find(opts.trajectory.idx==opts.idx(ii));
        if ~isempty(tix)
            scatter(ax(1),opts.trajectory.Xt(1,:),opts.trajectory.Zt(1,:),50,'w','LineWidth',opts.lw*1.5)
            plot(ax(1),opts.trajectory.Xt(1:tix,:),opts.trajectory.Zt(1:tix,:),'LineWidth',opts.lw*1.25)
            scatter(ax(1),opts.trajectory.Xt(tix,:),opts.trajectory.Zt(tix,:),50,'w','filled','MarkerEdgeColor','b','LineWidth',opts.lw*1.5)
        end
    end
    
    % Plot tracked objects?
    if ~isempty(opts.tracks)
        if isempty(opts.colorOrder)
            co = get(ax(1),'ColorOrder');
            co = co([1:2 4:7],:);
        else
            co = opts.colorOrder;
        end
%         count = 0;
        for tt=1:length(opts.tracks)
            [chk,trkIdx]=ismember(opts.idx(ii),opts.tracks(tt).tI);
            if chk
                cmask = false(length(z),length(x));
%                 count = count+1;
                cmask(opts.tracks(tt).clustIdx{trkIdx}) = true;
                [~,~,cpoly]=getROI(cmask,'maxRegions',1);
                
                ci = tt-floor(tt/size(co,1))*size(co,1);
                if rem(tt,size(co,1))==0; ci=size(co,1); end

                tp(tt) = plot3(ax(1),x(cpoly.X),z(cpoly.Y),TM+cpoly.Y*0,'Color',co(ci,:),'LineWidth',opts.lw*1.2);
            end
            
        end
    end
    
    if ~isempty(opts.roi)
        axis(ax(1),[opts.roi(3:4) opts.roi(1:2)])
    end
    
    if or( strcmp(opts.plotmode,'colorplot') , strcmp(opts.plotmode,'both') )
        axes(ax(2))
    end

    pause(opts.fdt)
end

if nargout>=1
    varargout{1} = ax;
end
if nargout>=2
    varargout{2} = ss;
end
if nargout>=3
    if ~isempty(opts.tracks) && exist('tp','var')
        varargout{3} = tp;
    else
        varargout{3} = plot(nan,nan);
    end
end
end