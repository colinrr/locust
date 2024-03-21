% PRE-PROCESS thermal imaages for plumeTracker where complex clouds etc
% cause problems for the tracking algorithm.
% ========================================================================
%  25A-0 step by step workflow - assumed to run from within workflowDriver
% clear all ; %close all;
% run project_25A0

   
%% STEP 1 - check 25A-4 polygons/masks to see if they work...
if get_foreground_mask
    verticalDisplace = 5; % Number of pixels to shift upwards to add a buffer pad ABOVE actual foreground
    poly_dc = [3 8]; % translation to apply to 25A4 vent clipping polygon


    polys_25A4.fileDir = fullfile(procDir,'25A4_preProcFiles/');
    polys_25A4.fgFile = fullfile(polys_25A4.fileDir,'foreground_mask_25A4.mat');
%     polys_25A4.fgPolyFile = fullfile(polys_25A4.fileDir,'foreground_poly_25A4.mat');
%     polys_25A4.histFile = fullfile(polys_25A4.fileDir,'Temp_time_hist_fgmask_25A4.mat');


    % Load image + 25A4 polygon
    load(heads)
    refImg = fullfile(matDir,T.File{num2str(ref_idx)});

    load(refImg)
    load(polys_25A4.fgFile)

    % Get mask
    fgMask = makeForegroundMask(refImg,thresh_fg);
    % Add buffer to cover foreground edge
    fgMask = [fgMask(verticalDisplace+1:end,:); ones(verticalDisplace,size(fgMask,2))];

    % Get mask from 25A-4 polygon for clipping out vent region
    ventClipMask = poly2mask(Poly.P(:,1)+ poly_dc(2),Poly.P(:,2)+ poly_dc(1), size(Frame,1), size(Frame,2));

    % Appy clip
    fgMask = and(fgMask, ~ventClipMask);
    
    save(fgFile,'fgMask','ventClipMask')

    
    if plot_flag
        figure
        subplot(1,2,1)
        imagesc(fgMask)
        title('Foreground Mask')
        subplot(1,2,2)
        
        histogram(Frame(:))
        hold on
        plot(thresh_fg*[1 1],ylim,'--','LineWidth',1.5)
        xlabel('Brightness Temperature (K)')
        ylabel('Counts')
        legend('Image','Foreground cutoff threshold')
        
        [~,fgMask,fgPoly] = getROI(fgMask,'maxRegions', 1);

        playFrames(matDir,heads,[],1,[],[],[230 330])
        hold on
        plot(fgPoly.X,fgPoly.Y,'LineWidth',1.5) % Mountain foreground poly
        plot(Poly.P(:,1) + poly_dc(2),Poly.P(:,2) + poly_dc(1),'LineWidth',1.5) % Rough non-plume exclusion poly
    end
end
%% STEP 2: Initial Temp-vs-time histograms

if view_histograms
    load(fgFile)
    load(heads)

    % Raw (no foreground mask)
    % H0 = getTimeHistogram(ptInputDir,ptInputHeads,[],[],[nullVal:330]);
    % figure
    Ht = getTimeHistogram(matDir,heads,[],fgMask,[nullVal:330]);

    % Plot Time Histograms
    if plot_flag
        figure
        % plotTimeHistogram(Ht) % Without threshold filter track
        plotTimeHistogram(Ht,ITW,T.Time) % With threshold filter track
    end
end

%% STEP 3: Apply masks and/or smooth heaviside filter - test first
if apply_filters_test
    load(fgFile)
    maskThermal(matDir,matHeads,procIdxTest,procDir,ITW,fgMask,nullVal)
    
    if plot_flag
        load(heads)
        % Compare histograms and frames of before and after filters
        F_pre = load(fullfile(matDir,T.File{num2str(ref_idx)}));
        F_post = load(fullfile(procDir,T.File{num2str(ref_idx)}));
        
        figure('Position',[100 100 1000 500])
        subplot(1,3,1)
        imagesc(F_pre.Frame)
        title('Raw image')
        cb = colorbar;
        cb.Label.String = 'Brightness T. (K)';
        subplot(1,3,2)
        imagesc(F_post.Frame)
        title('Filtered image')
        cb = colorbar;
        cb.Label.String = 'Brightness T. (K)';
        subplot(1,3,3)
        histogram(F_pre.Frame(:),100)
        hold on
        histogram(F_post.Frame(:),100)
        xlabel('Brightness T. (K)')
        ylabel('Counts')
        legend('Raw','Filtered')
        
        Ht = getTimeHistogram(procDir,heads,procIdxTest,fgMask,[nullVal:330]);
        figure
        plotTimeHistogram(Ht,ITW,T.Time(procIdxTest)') % With threshold filter track

        
    end
end


%% STEP 4: test plume_track
if flag_plumetrack_test
    %                                                             oDir,      ref, deb, fin, dt
    [content,refData] = mainTrackPlume(procDir,heads,plumeTrackDir, procIdxTest(1), procIdxTest, procIdxTest(end), dN, fixpix, PTplotFlags,delT);
    
    %.... fixed nullVal back to 210 as for 25A4
end

%% STEP 5: apply filter to all frames
if apply_filters
    load(fgFile)
    maskThermal(matDir,matHeads,[],procDir,ITW,fgMask,nullVal)
end
