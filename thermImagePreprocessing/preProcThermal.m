% PRE-PROCESS thermal imaages for plumeTracker where complex clouds etc
% cause problems for the tracking algorithm.
% ITERATIVE WORKFLOW:
%   1) Select REFERENCE FRAME (ref_idx) and FOREGROUND TEMPERATURE
%   THRESHOLD (thresh_fg) to create a foreground mask. Enter thermal image
%   SATURATION VALUE (satVal).
create_fg_mask = true;
%   2) Set initial temperature threshold (T0) and step width (dT) for
%   smooth heaviside filter. Plot histograms and frames for ref_idx, and
%   start, mid, end points of eruption
design_heaviside = false;
%   3) View tiled histograms to show effect of filter. Estimate values for
%   peakRange (temperature ranges containing the actual plume/object of
%   interest) in the frames.
%   4) Run initial filter, save output files, run plumeTracker.
%   5) Review output using pTfile and frame tile & hist tile plots
%   6) Use ITW to create a variable filter if needed.
%   7) Iterate 2-6 as needed.

% clear all; close all
% Pre-process thermal frames 
  % Things to try:  thresholding...actually...meh maybe not
  %                 gradient'

  clear all; close all
% Get input folders and parameters for which event?
% EVENT = '24A';
% EVENT = '24B';
EVENT = '25A4';
% EVENT = '25B';

run preProcInput_25A0
% project_25A4

%%
% function [something] = preProcThermal(dirS,parS,PS,flags,dd)
%   INPUT:    All inputs are structures
%           dirS = directories: matDir, procDir, heads, PTfile
%           parS = parameters: 
%                         thresh_fg
%                         satVal
%                         nullVal
%                         T0
%                         dT
%                         peakRange
%                         ITW
%                         interp_meth
%                         bins
%                         ref_idx
%                         Idx
%
%           PS   = polygon structure
%                         design_polys
%                         PolyFile
%                         show_polys
%                         npolys
%                         IdxPoly
%                         DesPol_Idx
%
%           flags:        design_flag % Use this when designing the filters from a reference image - allows plotting in separate figures
%                         plot_flag  % Use this when designing or vetting, turn off for mass production...or turn on and save?...
%                           imstat   % Plot filter, histogram, and images steps for each file
%                           plot_idx % Select indices to plot imstat\
%                           hist_tile % Tile histograms for each file
%                           fram_tile % Tile RAW frames...with outlines if provided
%                           plot_idx2 % Select indices for which to plot hist and frame tiles
%                         write_flag  % Save output, overwrite existing frames in procDir
%           
%           dd: display params
%                         rows  % N rows and columns for tiled plotted
%                         cols
%                         dd_filt % figure position vectors
%                         dd_hist
%                         dd_mask
%                         dd_rawi
%                         dd_mski
%                         dd_scli
%% DO THE THING ================================================

if ~isempty(ITW)
    fprintf('Using variable thresholding...\n')
%     Do some stuff here..
    T0v = interp1(ITW(1,:),ITW(2,:),Idx,interp_meth);
    dTv = interp1(ITW(1,:),ITW(3,:),Idx,interp_meth);
    Kv = 2./dTv;
end

% Fig params
pdx = 0.05;
pdy = 0.08;
% design dimensions

load(heads)
% n= numel(thresh_fg);

hs_sm = (@(T,T0,k) 1./(1 + exp(-2*k*(T-T0))));
Ttest = linspace(nullVal,satVal,201);
Te    = bins(1:end-1) + diff(bins);
C1 = zeros(numel(Te),numel(Idx)); % Matrix of histcounts
C2 = C1;
C3 = C1;

% ----- 1) Get foreground mask from reference image -------
load(fullfile(matDir,T.File{num2str(ref_idx)}))
Fgmask         = Frame>thresh_fg;

if ~isempty(pTfile); ptT=load(pTfile); Tout=ptT.T; end


%% ================= POLYGON DESIGN ====================
if design_polys
    np = numel(npolys);
    pf = figure('position',[50 200 1800 1000]);
    
%     for nn=1:numel(DesPol_Idx)
%         load(fullfile(matDir,T.File{num2str(DesPol_Idx(nn))}))
%         pax(nn) = tightSubplot(1,numel(DesPol_Idx),nn,0);
%         if ismember(num2str(DesPol_Idx(nn)),Tout.Properties.RowNames)
%             imagesc(Frame + Frame.*Tout{num2str(DesPol_Idx(nn)),'Outline'}{1});
%         else
%             imagesc(Frame)
%         end
%         axis off
%         hold on
%         caxis([190 350])
% %         axis equal
%         daspect([1 1 1])
% 
%     end
%     colormap(thermgray(150))
    
    for nn=1:np
        despol_idx = DesPol_Idx{nn};
        for mm=1:numel(despol_idx)
                load(fullfile(matDir,T.File{num2str(despol_idx(mm))}))      
                pax(mm) = tightSubplot(2,2,mm,0,0);
            if ismember(num2str(despol_idx(mm)),Tout.Properties.RowNames)
                imagesc(Frame + Frame.*Tout{num2str(despol_idx(mm)),'Outline'}{1});
            else
                imagesc(Frame)
            end
            axis off
            hold on
            caxis([190 350])
    %         axis equal
            daspect([1 1 1])

        end
        colormap(thermgray(150))           
            
        pf2=figure('position',[50 200 1500 1000]);
        load(fullfile(matDir,T.File{num2str(despol_idx(1))}))
        if ismember(num2str(despol_idx(1)),Tout.Properties.RowNames)
            imagesc(Frame + Frame.*Tout{num2str(despol_idx(1)),'Outline'}{1});
        else
            imagesc(Frame)
        end
        axis off
        hold on
        colormap(thermgray(150))
        caxis([190 350])
%         linkaxes([pax gca],'xy')
        axis equal        

        polygon = getPoly([pax gca],size(Frame));
%         if ~npolys(nn)
%             polygon = ~polygon;
%         end
        Polys(nn).P          = polygon;
        Polys(nn).Type       = npolys(nn);
        Polys(nn).Design_Idx = despol_idx;
        Polys(nn).Idx        = IdxPoly{nn};
        Polys(nn).Mask       = poly2mask(polygon(:,1),polygon(:,2),size(Frame,1),size(Frame,2));
        
        if ~Polys(nn).Type
            Polys(nn).Mask   = ~Polys(nn).Mask;
        end
%         Polys(nn).Outline    = sparse(edge(polygon,'sobel')); 
        

        clf(pf)
        clf(pf2)
    end
    
    % Convert poly to mask and add to a cell/struct array
    
    save(PolyFile,'Polys')
    close(pf)
    close(pf2)
end


%% ================ Apply filters and plot ================
% Plotting of tiled histograms
if hist_tile;  hf = figure('position',[50 50 1500 1000]); hcount=0; end
if fram_tile; ff = figure('position',[50 1050 1500 1000]); fcount=0; end
if ~isempty(peakRange); peakT = Idx*0; peakN = Idx*0; tf = figure('position',[50 50 1200 600]); end
if and(show_polys,~design_polys); load(PolyFile); end

% MANUAL EDIT
% Polys.Idx = 928:1576;

if write_flag
    fprintf('Writing processes frames to:\n\t%s\n',procDir)
end
for ii=1:numel(Idx)
    idx = Idx(ii);
    
    load(fullfile(matDir,T.File{num2str(idx)}))
    
    % ------ 2) Apply foreground mask --------
    Frame2        = Frame;
    Frame2(Fgmask) = nullVal;



    % ------ 3) Apply smooth heaviside filter ------
    if ~isempty(ITW) % Use non-constant filter?
        T0 = T0v(ii);
        dT = dTv(ii);
        K  = Kv(ii);
    else

        K = kw;
    end
    
    Frame3  = hs_sm(Frame2,T0,K) .* (Frame - nullVal) + nullVal;
 
    N  = histcounts(Frame,bins);
    N2 = histcounts(Frame2,bins);
    N3 = histcounts(Frame3,bins);
    
    C1(:,ii) = N;
    C2(:,ii) = N2;
    C3(:,ii) = N3;
   
    % Find temperature of for peak count in plume
    if ~isempty(peakRange)
        [~,pI]  = closest(idx,peakRange(1,:),'prev');
        rgI     = and(Te>=peakRange(2,pI),Te<=peakRange(3,pI));
        Tr      = Te(rgI);
        [peakN(ii),pkI] = max(N2(rgI));
        peakT(ii)  = Tr(pkI);
    end
    
% -------- Not very useful gradient processing...
%     [Gmag,Gdir] = imgradient(Frame);
    
%     F = Frame./T_scale;
%     G = Gmag./dT_scale;
    
    %Blur
%     Fb = imgaussfilt(F,2);

    % Composite frame
%     Frame = F+G*rel_scale;


%         figure('position',[0 800 1600 500])

 %% ================ PLOTTING ==========================
    if ii==1
        fn = sprintf('Reference image: %i',idx);
    else
        fn = sprintf('Frame index: %i',idx);
    end
    
    if plot_flag
        Tm = max([max(Frame(:)) max(Frame2(:)) max(Frame3(:))]);

        if imstat && ismember(idx,plot_idx)
            % Plot filter function
            if design_flag
                figure('position',dd_filt,'Name',fn)
            else
                figure('position',[50 80 1500 1000],'Name',fn)
                tightSubplot(2,3,1,pdx,pdy);
            end
            plot(Ttest,hs_sm(Ttest,T0,K),'linewidth',1.7) %*(satVal-Tmin)+Tmin)
            hold on
            plot((T0-dT)*[1 1],[0 1],'--k')
            plot((T0+dT)*[1 1],[0 1],'--k')
            title(sprintf('Filter: T0=%.1f, dt=%i, k=%.2f',T0,dT,K))
            xlabel('Kelvin')
            axis tight

            % Plot histograms
            if design_flag
                figure('position',dd_hist)
            else
                tightSubplot(2,3,2,pdx,pdy);
            end
            fill(T0 + dT*[-1 1 1 -1], [10^0 10^0 max([N N2 N3]) max([N N2 N3])],[1 0.95 0.85],'EdgeColor','None')
            hold on
            plot(Te,N,'.-b','linewidth',1.4)
    %         hold on
            plot(Te,N2,'.-g')
            plot(Te,N3,'.-r')
            % plot([1 numel(Fnew)],T0+dT,'--k')
            set(gca,'YScale','log')
            title('Histogram counts')
            xlabel('Kelvin')
            axis tight
            legend({'Scale range','Raw','Masked','Scaled'})

            % Foreground mask
            if design_flag
                figure('position',dd_mask)
            else        
                tightSubplot(2,3,3,pdx,pdy);
            end
            imagesc(Fgmask)
            axis off
    %         caxis([Tmin Tm])
            axis equal
            title('Mask')

            % Raw image
            if design_flag
                figure('position',dd_rawi)
            else
                tightSubplot(2,3,4,0,pdy);
            end
            imagesc(Frame)
            axis off
            caxis([nullVal Tm])
            axis equal
            title('Raw image')

            % Masked image
            if design_flag
                figure('position',dd_mski)
            else
                tightSubplot(2,3,5,0,pdy);
            end
            if ~isempty(pTfile) && ismember(num2str(idx),Tout.Properties.RowNames)
                imagesc(Frame2 + Frame2.*Tout{num2str(idx),'Outline'}{1});
            else
                imagesc(Frame2)
            end
            axis off
            caxis([nullVal Tm])
            axis equal
            title('Masked')

            % Masked and scaled image
            if design_flag
                figure('position',dd_scli)
            else
                tightSubplot(2,3,6,0,pdy);            plot(thresh_fg*[1 1],[1 10^5],'--k')

            end
            if ~isempty(pTfile) && ismember(num2str(idx),Tout.Properties.RowNames)
                imagesc(Frame3 + Frame3.*Tout{num2str(idx),'Outline'}{1});
            else
                imagesc(Frame3)
            end
            axis off
            caxis([nullVal Tm])
            axis equal
            title('Scaled and masked')
            colormap(CubeHelix)
        end

        % Tiled histograms
        if hist_tile && ismember(idx,plot_idx2)
            hcount = hcount+1;
            figure(hf)
            ax=tightSubplot(rows,cols,hcount,0,0);
            fill(T0 + dT*[-1 1 1 -1], [10^0 10^0 max([N N2 N3]) max([N N2 N3])],[1 0.95 0.85],'EdgeColor','None')
            hold on
            ax.ColorOrderIndex=1;
%             set(gca,'ColorOrder','factory')
            plot(Te,N2,'linewidth',1.6)
            plot(Te,N3)
            plot(thresh_fg*[1 1],[1 10^5],'--k')
            set(gca,'YScale','log')
            if ~ismember(ii,(0:rows-1)*cols+1)
                set(gca,'YTickLabel',[])
            end
            if ii<=((rows-1)*cols)
                set(gca,'XTickLabel',[])
            end
            axis tight
            xlim([min(bins) max(bins)])
            text(0.7,0.85,sprintf('%i',idx),'Units','normalized')
            
            if ~isempty(peakRange)
                scatter(peakT(ii),peakN(ii),10,'ok','filled')
            end
        end
        
        % Tiled frames
        if fram_tile && ismember(idx,plot_idx2)
            fcount = fcount +1;
            figure(ff)
            ax=tightSubplot(rows,cols,fcount,0,0);
%             tightSubplot(rows,cols,ii,0,0);
            if ~isempty(pTfile) && ismember(num2str(idx),Tout.Properties.RowNames)
                imagesc(Frame + Frame.*Tout{num2str(idx),'Outline'}{1});
            else
                imagesc(Frame)
            end
            caxis([nullVal 370])
            set(gca,'XTickLabel',[],'YTickLabel',[])
            colormap(thermgray(100))
            text(0.7,0.85,sprintf('%i',idx),'Units','normalized','Color',[0.9 0.9 0.9])
            if show_polys
                for pp=1:length(Polys)
                    if ismember(idx,Polys(pp).Idx)
                        hold on
                        if Polys(pp).Type
                            ecol = [0.9 0.9 0.9];
                            fcol = ecol;
                            alph = 0.3;
                        else
                            ecol = [0 0 0];
                            fcol = [0.2 0.2 0.2];
                            alph = 0.4;
                        end
                        fill(Polys(pp).P(:,1),Polys(pp).P(:,2),fcol,'edgecolor',ecol,...
                            'facealpha',alph);
                    end
                end
            end
            if ~isempty(pTfile) && ismember(num2str(idx),Tout.Properties.RowNames)
                xx = sum(Tout{num2str(idx),'Outline'}{1},1);
                yy = sum(Tout{num2str(idx),'Outline'}{1},2);
                xlim([find(xx,1,'first') find(xx,1,'last')])
                ylim([find(yy,1,'first') find(yy,1,'last')])
            end

        end
    end


% ----- WRITING ------
    if write_flag
    % Save gradient file
        Frame = Frame3;
        oname = fullfile(procDir,T.File{num2str(idx)});
        fprintf('\t%i:\t%s\n',idx,T.File{num2str(idx)})
        save(oname,'Frame','File_DateTime')

    end
end

if write_flag
        disp('Writing pre-processing paramter file.')
        preProc.satVal      = satVal;
        preProc.nullVal     = nullVal;
        preProc.thresh_fg   = thresh_fg;
        preProc.ref_idx     = ref_idx;
        preProc.Idx         = Idx;
        preProc.T0          = T0;
        preProc.dT          = dT;
        preProc.peakRange   = peakRange;
        preProc.ITW         = ITW;
        preProc.interp_meth = interp_meth;
        preProc.bins        = bins;

        if ~isempty(ITW)
            preProc.T0v = T0v;
            preProc.dTv = dTv;
            preProc.Kv  = Kv;
        end
        save(fullfile(procDir,'preProcHeads.mat'),'preProc');
end





if ~isempty(peakRange) && plot_flag
    figure(tf)
    legs={};
    labs={};
    if size(C1,2)>=20 && plot_flag
        axa=tightSubplot(1,3,1,0);
        pcolor(Idx,Te,log10(C1));
        shading flat
        title('Raw images')
        colormap(gray)

        axb=tightSubplot(1,3,2,0);
        pcolor(Idx,Te,log10(C2));
        shading flat
        colormap(gray)
        title('Masked images')

        axc=tightSubplot(1,3,3,0);
        pcolor(Idx,Te,log10(C3));
        shading flat
        colormap(gray)
        title('Scaled images')

        linkaxes([axa axb axc], 'y')
        lax = axa;

        axes(axb)
        hold on
    
    end
    
    aa=plot(Idx,peakT,'.-');
    hold on
    legs = [legs aa];
    labs = [labs 'Apparent median plume temp'];
end
if ~isempty(ITW) && plot_flag
    plot(ITW(1,:),ITW(2,:),'o')
    hold on
    ab=plot(Idx,T0v);
    ac=plotLineError(Idx,T0v,dTv,ab.Color,0.25);
    legs = [legs ab ac];
    labs = [labs 'Heaviside temp threshold' 'Heaviside step range'];
end
if and(or(~isempty(peakRange),~isempty(ITW)),plot_flag)
    if size(C1,2)>=30; axes(lax); end
    xlabel('Frame Index')
    ylabel('Kelvin')
    legend(legs,labs)
end

%%
function Pout = getPoly(h,sz)
% Function for building polygons
% h = axes handles
% sz = size(Frame)
b = 1;
X = [];
Y = [];
P = []; % plot handles
for hh=1:numel(h)
    P = [P plot(1,1)]; % dummy handles
end

disp('Design polygon mask.')
disp('Commands:')
disp('   left click: add point')
disp('   right click: delete point')
disp('   z: zoom 2x towards cursor')
disp('   x: zoom out 2x away from cursor')
disp('   y: quit and save poly')
disp('   q: quit and cancel poly')


while and(b~=121,b~=113)
    axes(h(end))
    [x,y,b] = ginput(1);
    axl = axis;

    % THE CLICKING BIT
    % Commands: 
    %   left click: add point
    %   right click: delete point
    %   z: zoom 2x towards cursor
    %   x: zoom out 2x away from cursor
    %   y: quit and save poly
    %   q: quit and cancel poly

    if b==1 % left click
        X = [X x];
        Y = [Y y];
    elseif b==3 % right click
        X = X(1:end-1);
        Y = Y(1:end-1);
    elseif b==122 % z
        dx = diff(axl(1:2));
        dy = diff(axl(3:4));
%         x0 = mean(axl(1:2));
%         y0 = mean(axl(3:4));      
        axl = [x-dx/4 x+dx/4 y-dy/4 y+dy/4]; 
    elseif b==120 % x
        axl = axis;
        dx = diff(axl(1:2));
        dy = diff(axl(3:4));
        x0 = mean(axl(1:2));
        y0 = mean(axl(3:4));
        axl = [x-dx x+dx y-dy y+dy];         
    elseif b==121 % y
        save_flag = true;
    elseif b==113 % q
        save_flag = false;
    end
    
    for hh = 1:numel(h)
            axes(h(hh))
        if or(b==122,b==120)
           axis(axl)
        end
        if or(b==1,b==3)
            delete(P(hh))
            P(hh) = plot(X,Y,'o-w');
        end
    end
end

if save_flag
    Pout = [X' Y']; %poly2mask(X,Y,sz(1),sz(2));
else
    disp('Selection cancelled, polygon discarded.')
    Pout = [];
end

end