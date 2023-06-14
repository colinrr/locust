function [f1,f2,tax,avax] = plotScalingResultsMS(tk,dat,events,avgErrMode,plotByTime,plotR0,oPath)

%close all
% INPUT:
% dat
% tk?
% sI
% events = {'25A4','24A','25B'};
% xmode? = event number or time?

overlaysOff = true;

if nargin<4
    avgErrMode = 2;
end
if nargin<5 || isempty(plotByTime)
    plotByTime = true;
end
if nargin<6 || isempty(plotR0)
    plotR0 = false;
end

dataDir   = '~/Kahuna/data/sabancaya_5_2018/';  % Calculon
% avgImgFile = fullfile(dataDir,'pulseTrack_analysis/AverageImageScalingResults.mat');
% avgImgFile = fullfile(dataDir,'pulseTrack_analysis/AverageImageScalingResults_V2.mat');
avgImgFile = fullfile(dataDir,'pulseTrack_analysis/AverageImageScalingResults_V3.mat');
thermSources = fullfile(dataDir,'pulseTrack_analysis/allSources_2020-10-29_24A_25A4_25B.mat');
avI = [1 2 3]; % Order of events?

load(avgImgFile)
load(thermSources)

nd = length(dat); % Number of runs
ne = length(events);


% for ti = 1:nt
%     tX(ti) = avX(evNo(ti)) + lasttx(evNo(ti));
%     lasttx(evNo(ti)) = lasttx(evNo(ti))+1;
% end
% tXl = string(sIc.trackSet);
% boundX = avX(1:2) + trksPerEvent(1:2) +0.5;

%% Plotting params and settings
% Plot track scaling results
% --> drdz, r0, z0, B
if plotR0
    nc = 5;
else
    nc = 4;
end

nr = ne;
% nc = 4; 
dx = 0.01;
dy = 0.013;
ppads = [0.1 0.04 0.08 0.08];
cbY = 0.93; % Color bar x position
ySz = [];
% ySz = [1 1 1 1.5];
fpos  = [50 50 900 900];
lpos  = [0.1 cbY 0.1 0.025];
odims = [18 16];
fs = 14;
lfs = 13;
lw = 1.1;
lw2 = 1;

% For print
fs = 9;
lfs = 8;

avSym = 'o';
% tkAvSym = 's';
runSym = {'s','d','^','v','p','h'};
dxs = 0.1; % Plot shift for subsequent data points

tkCol = [0.3 0.3 0.3];
tkCol2 = [0.4 0.4 0.4];
alph = 0.25;
prctileI = 4;


% cmap = redblue; %'winter';
% cmap = parula;
cn = winter;
cn = cn(1,:);
% c1 = [1 1 1];
c1 = [1 0 0];
ncol = 150;
cmap = flipud([linspace(c1(1),cn(1),ncol)' linspace(c1(2),cn(2),ncol)' linspace(c1(3),cn(3),ncol)']);


% Positions for t-avg and track-average values
if plotByTime
    avX = -10;
    tkAvX = -5;
else
    avX = 1;
    tkAvX = 2;
end

% Y Limits for time plots
% customyl = [avX 300
%             avX 160
%             avX 70];

dt25 = [0 0 16.646]; % Time delay from 25B event start to video start
customyl = [0 300
            0 160
            -2 70+dt25(3)];
customxl = [-700 375
            0 0.58
            -4.6 0.1];
% customxl = [-790 375
%             0 0.58
%             -5.5 0.2];
BTick = [-5:1:0];
drTick = [0:0.1:0.5];
%%

% Fit error measure limits
    rrmseAll = []; %zeros(nd*nt,1);
    nrmseAll = []; %zeros(nd*nt,1);
    for di=1:nd
            rrmseAll = [rrmseAll; dat(di).fitArrays.r_rmse];
            nrmseAll = [nrmseAll; dat(di).fitArrays.B_nrmse];
    end
    rcl = [nanmin(rrmseAll) nanmax(rrmseAll)];
    bcl = [nanmin(nrmseAll) nanmax(nrmseAll)];
    if numel(unique(bcl))==1; bcl = [0 bcl(1)]; end


% Axes positions
if nc==4
    ai0=0;
    xSz = [0.8 1 1 1.2];
elseif nc==5
    ai0=1;
    xSz = [0.8 1 1 1 1.2];
else
    error('What did you do?')
end
    
% Order of average images
%% Results plot

% Setup figure
f1=figure('position',fpos);
taxxl = zeros(ne,2);
for ei = 1:ne
    tax(ei)    = tightSubplot(nr,nc,1+nc*(ei-1),dx,dy,ppads,xSz,ySz);
    sei = find(ismember({S.event},events{ei}));
    plot(S(sei).T0.prctile(prctileI,:)./max(S(sei).T0.prctile(5,:)),S(sei).T0.t+dt25(ei),'k','LineWidth',lw)
    hold on
    plot(S(sei).T0.prctile(5,:)./max(S(sei).T0.prctile(5,:)),S(sei).T0.t+dt25(ei),'Color',[0.5 0.5 0.5],'LineWidth',lw)
%     taxxl(ei,:) = [0 max(S(sei).T0.prctile(5,:))];
    taxxl(ei,:) = [0 1];
    axis tight
    xlim(taxxl(ei,:))
    ylim(customyl(ei,:))
    grid on
    cot = get(tax(ei),'ColorOrder');

    for ci = 1:(nc-1)
        avax(ei,ci) = tightSubplot(nr,nc,ci+1+nc*(ei-1),dx,dy,ppads,xSz,ySz);
%         avax(ei,2) = tightSubplot(nr,nc,3+nc*(ei-1),dx,dy,ppads,xSz,ySz);
%         avax(ei,3) = tightSubplot(nr,nc,4+nc*(ei-1),dx,dy,ppads,xSz,ySz);
%     avax(ei,4) = tightSubplot(nr,nc,5+nc*(ei-1),dx,dy,ppads,xSz,ySz);
    end

    hold(avax(ei,:),'on')

    
% Plot averaged image results
    if plotByTime
%         if nc==5
%             axes(avax(ei,1))
%             plotLineError(AvgImgT.r0_combo(avI(ei),:).*[1;1],TavgS(avI(ei)).tSpan',tkCol2,alph)
%         end
        axes(avax(ei,ai0+1))
        if ei==3 % Adjusting 25B plot to account for the early start time of track 1
            tspan = [dt25(ei) TavgS(avI(ei)).tSpan(2)]';
        else
            tspan = TavgS(avI(ei)).tSpan'; 
        end
        plotLineError(AvgImgT.z0_combo(avI(ei),:).*[1;1],tspan,tkCol2,alph)
        axes(avax(ei,ai0+2))
        plotLineError(AvgImgT.drdz_combo(avI(ei),:).*[1;1],tspan,tkCol2,alph)
        axes(avax(ei,ai0+3))
        plotLineError(AvgImgT.B_combo(avI(ei),:).*[1;1],tspan,tkCol2,alph)        
    else    
        if nc==5
            axes(avax(ei,1))
            [sh,eh] = plotScaleError(AvgImgT.r0_combo(avI(ei),:),avX,0,[],[],avSym);
        end
        axes(avax(ei,ai0+1))
        [sh,eh] = plotScaleError(AvgImgT.z0_combo(avI(ei),:),avX,0,[],[]);
        axes(avax(ei,ai0+2))
        [sh,eh] = plotScaleError(AvgImgT.drdz_combo(avI(ei),:),avX,0,[],[]);
        axes(avax(ei,ai0+3))
        [sh,eh] = plotScaleError(AvgImgT.B_combo(avI(ei),:),avX,0,[],[]);
    end
end
hold(tax,'on');
set(tax,'FontSize',fs,'YDir','reverse')
xlabel(tax(ne),'$\displaystyle{\Delta T/\Delta T_{max}}$ (\%$_{ile}$)','Interpreter','Latex')
linkaxes(tax,'x')
set(tax(1:2),'XTickLabel',[])
%% Plot source time-series

%% Plot track results

xl = zeros(ne,2,nd);
tkYtick = cell(ne,1);
tkYtickL = cell(ne,1);
% yl = zeros(nd,2);
for di = 1:nd % Loop over runs
    
    nt = length(dat(di).in.trackSet);
    evNo = zeros(nt,1);
    tSpan = zeros(nt,2);
    for ti = 1:nt
        if isfield(dat(di).fitArrays,'t1_srcWin')
            tSpan(ti,1) = dat(di).fitArrays.t1_srcWin(ti);
        else
            tSpan(ti,1) = tk(dat(di).in.trackSet(ti)).t(1);
        end
        tSpan(ti,2) = tk(dat(di).in.trackSet(ti)).t(end);
        [~,evNo(ti)] = ismember(dat(di).trackFits(ti).eventName,events);
%         [~,evNo(ti)] = ismember(tk(dat(di).in.trackSet(ti)).event,events);
    end
    trksPerEvent = zeros(1,ne);
    for ei = 1:ne
        trksPerEvent(ei) = (sum(evNo==ei));
    end
    
    % Loop over events
    for ei = 1:ne
        tkIdx = (evNo==ei); % Which tracks belong to this event?
        Rqc = [dat(di).trackFits(tkIdx).Rqc];
        Tqc = [dat(di).trackFits(tkIdx).Tqc];
        Tqc(isnan(dat(di).fitArrays.B_nrmse(tkIdx,:)))=0;
        
        % Plot by event number
%         tkX = tkAvX + [tk(dat(di).in.trackSet(tkIdx)).eventTrack]; % X positions for track subset
        
        % Plot by time
        if isfield(dat(di).fitArrays,'t1_srcWin')
            tkX = dat(di).fitArrays.t1_srcWin(tkIdx);
        else
            tkX = tSpan(tkIdx,1);
        end
        if ei==3 % Apply a manual adjustment to 25B track 1 to reflect its true start time
            tkX(1) = 0;
            tkX(2:end) = tkX(2:end) + dt25(ei);
        end
        tp{ei} = plot(tax(ei),repmat(taxxl(ei,:)',[1 length(tkX)]),(tkX.*[1 1])','--','LineWidth',lw,'Color',cot(1,:));
        
%         if nc==5
%             % r0
%             axes(avax(ei,1))
%             if any(tkIdx)
%                 % Mean track results
%                 meanY = yErrMu(dat(di).fitArrays.r0(tkIdx,:),dat(di).fitArrays.r_rmse(tkIdx),avgErrMode);
% %                 [sh,eh] = plotScaleError(meanY,tkAvX+ dxs.*(di-1),mean(dat(di).fitArrays.r_rmse(tkIdx,:)),[],[],runSym{di});
%                 %Indiviual track results
%                 [sh,eh] = plotScaleError(dat(di).fitArrays.r0(tkIdx,:),tkX + dxs.*(di-1),dat(di).fitArrays.r_rmse(tkIdx,:),Rqc,[],runSym{di},tkCol);
%             end
%             if ei==ne
%     %             title(sprintf('Event %i',ei))
%                 caxis(rcl)
%             end
%         end
        
        % z0
         axes(avax(ei,ai0+1))
        if any(tkIdx)
            % Plot average of tracks
            meanY = yErrMu(dat(di).fitArrays.z0(tkIdx,:),dat(di).fitArrays.r_rmse(tkIdx),avgErrMode);
            plotLineError(meanY.*[1;1],tkX([1 end]),cot(di,:),alph,true);
%             [sh,eh] = plotScaleError(meanY,tkAvX+ dxs.*(di-1),mean(dat(di).fitArrays.r_rmse(tkIdx,:)),[],[],runSym{di});
            % Plot track results
%             plot(repmat(customxl(1,:)',[1 length(tkX)]),(tkX.*[1 1])','--','LineWidth',lw2,'Color',tkCol2);
            [sh,eh] = plotScaleError(dat(di).fitArrays.z0(tkIdx,:),tkX + dxs.*(di-1),dat(di).fitArrays.r_rmse(tkIdx,:),Rqc,[],runSym{di},tkCol);
        end
        if ei==ne
            caxis(rcl)    
        end
        
        % drdz
        axes(avax(ei,ai0+2))
        if any(tkIdx)
            % Plot average of tracks
            meanY = yErrMu(dat(di).fitArrays.drdz(tkIdx,:),dat(di).fitArrays.r_rmse(tkIdx),avgErrMode);
            plotLineError(meanY.*[1;1],tkX([1 end]),cot(di,:),alph,true);
%             [sh,eh] = plotScaleError(meanY,tkAvX+ dxs.*(di-1),mean(dat(di).fitArrays.r_rmse(tkIdx,:)),[],[],runSym{di});
            % Plot track results
%             plot(repmat(customxl(2,:)',[1 length(tkX)]),(tkX.*[1 1])','--','LineWidth',lw2,'Color',tkCol2);
            [sh,eh] = plotScaleError(dat(di).fitArrays.drdz(tkIdx,:),tkX + dxs.*(di-1),dat(di).fitArrays.r_rmse(tkIdx,:),Rqc,[],runSym{di},tkCol);
        end
        if ei==ne
            caxis(rcl)
        end
    
        % B
        axes(avax(ei,ai0+3))
        if any(tkIdx)
            % Plot average of tracks
            meanY = yErrMu(dat(di).fitArrays.Bpl(tkIdx,:),dat(di).fitArrays.B_nrmse(tkIdx),avgErrMode);
            plotLineError(meanY.*[1;1],tkX([1 end]),cot(di,:),alph,true);
%             [sh,eh] = plotScaleError(meanY,tkAvX+ dxs.*(di-1),nanmean(dat(di).fitArrays.B_nrmse(tkIdx,:)),[],[],runSym{di});
            % Plot track results
%             plot(repmat(customxl(3,:)',[1 length(tkX)]),(tkX.*[1 1])','--','LineWidth',lw2,'Color',tkCol2);
            [sh,eh] = plotScaleError(dat(di).fitArrays.Bpl(tkIdx,:),tkX + dxs.*(di-1),dat(di).fitArrays.B_nrmse(tkIdx,:),Tqc,[],runSym{di},tkCol);
        end
        if ei==ne
            caxis(bcl)
        end
        
        
        % Ticks & labels
%         tkYtick{ei} = union(tkYtick{ei},[avX; tkAvX; tkX]);
        tkYtick{ei} = union(tkYtick{ei},[tkX]);
%         tkXtick{ei} = union(tkXtick{ei},dat(di).in.trackSet(tkIdx)); % Overall track No.
        tkYtickL{ei} = union(tkYtickL{ei},[tk(dat(di).in.trackSet(tkIdx)).eventTrack]); % Within-event track No.
%         xl(ei,:,di) = 
    end
end

for ei=1:ne
    set(avax(ei,1:(nc-2)),'YTickLabel',[])
    set(avax(ei,1:(nc-1)),'YTick',tkYtick{ei},'YLim',[0.5 tkYtick{ei}(end)+0.5],'FontSize',fs)
%     ytl = [{"Track" "Average"}; string(tkYtickL{ei})];
    ytl = string(tkYtickL{ei});
    set(avax(ei,nc-1),'YTickLabel',ytl) %,'YAxisLocation','right')
%     set(avax(ei,nc-1),'YTickLabel',[])
    for ai = 1:(nc-1)
        colormap(avax(ei,ai),cmap)
    end
    if plotByTime
        linkaxes([tax(ei) avax(ei,:)],'y')
        yl = customyl(ei,:);
    else
        yl = [0.5 tkYtick{ei}(end)+0.5];
    end
    if ~overlaysOff
        plot(avax(ei,nc-1),[-5/3 -5/3],yl,'--k','LineWidth',lw)
        plot(avax(ei,nc-1),[-3 -3],yl,':k','LineWidth',lw)
    end
    ylabel(tax(ei),{sprintf('\\bf{Event %i}',ei),'$t$ (s)'},'Interpreter','Latex')
end
   
cbStrings = {'NRMSE','NRMSE','NRMSE'};
for ai = 1:(nc-1)
    cb=colorbar(avax(1,ai),'northoutside');
    cb.Label.String = cbStrings{ai};
    cbpos = cb.Position;
    cb.Position = [cbpos(1) cbY cbpos(3:4)];
    linkaxes(avax(:,ai),'x')
    xlim(avax(:,ai),customxl(ai,:))
end
set(avax(:,nc-2),'XTick',drTick)
set(avax(:,nc-1),'YAxisLocation','right','XTick',BTick)
set(avax,'YDir','reverse')
set(avax(1:(ne-1),:),'XTickLabel',[])
% xlim(avax(4,:),[-5 0.5])
if nc==5
    xlabel(avax(ne,1),'$r_0$ (m)','Interpreter','Latex')
end
xlabel(avax(ne,ai0+1),'$\displaystyle{z_0}$ (m)','Interpreter','Latex')
xlabel(avax(ne,ai0+2),'$\displaystyle{\frac{dR}{dz}}$','Interpreter','Latex')
xlabel(avax(ne,ai0+3),'$\displaystyle{B}$','Interpreter','Latex')

grid(avax(:),'on')
box(avax(:),'on')

% Dummy legend
dp(1) = plot(nan,nan,'Color',tkCol2,'LineWidth',lw);
dp(2) = plot(nan,nan,'Color',cot(1,:),'LineWidth',lw);
[dp(3),~] = plotScaleError([nan nan nan],nan,0,1,gca,'s',tkCol);
% Track
% Poor track?
if ~overlaysOff
    legend(avax(1),dp,{'Time-average','Track Average','Track Result'},'Position',lpos,'Interpreter','Latex')
end
    %% B vs z0
    clear sh
    xl = [-450 350];
    f2 = figure;
    sh(1)=plot(xl,[-5/3 -5/3],'--k','LineWidth',lw);
    hold on
    sh(2)=plot(xl,[-3 -3],':k','LineWidth',lw);
    ll = {'Pure Plume','Pure Thermal'};
    co = get(gca,'ColorOrder');
    ct = 0;
    for di= 1
        for ei = 1:ne
            ct = ct+1;
            tkIdx = (evNo==ei);
            Tqc = [dat(di).trackFits(tkIdx).Tqc];
            Rqc = [dat(di).trackFits(tkIdx).Rqc];
            [sh(ct)]=plotScaleError(dat(di).fitArrays.z0(tkIdx,:),dat(di).fitArrays.Bpl(tkIdx,:),dat(di).fitArrays.r_rmse(tkIdx,:),Rqc,gca,runSym{di},co(ei,:));
            ll{ct} = sprintf('Event %i',ei);
        end
    end
    legend(sh,ll,'Location','northwest')
    xlabel('$z_0$ (m)','Interpreter','Latex')
    ylabel('$B$','Interpreter','Latex')
    grid on
    axis tight
    colormap(cmap)
    caxis(rcl)
    cb = colorbar('location','northoutside');
    cb.Label.String = 'NRMSE';

end