function varargout = regFrameDiff(matDir,matHeads,regDir,regHeads,roi,refidx)
% regFrameDiff(matDir,matPar,regDir,regPar,roi,refidx)
% Simple script to calculate mean and max differences between indiviudal 
% .mat file frames. This is useful to detect sequences that have missing
% time chunks, or where frames move and need to be adjusted with pixel
% registration. Can also use to check frames after applying registration.
%
% Scanning through all frames in the "paramf" data table, this function 
% calculates the mean and max absolute value of the difference between
% frames, using either a fixed reference frame, or between two
% consecutive frames.
%    matDir = directory to .mat frame files.
% OPTIONAL:
%    matPar = path to parameter table containing frame meta-data
%    regDir = optional path to registered directory for direct comparison
%           of before/after registration. Leave empty [] to ignore.
%    regPar = parameter table for registered frames in regDir. Leave empty
%           to ignore.
%    roi    = pixels limits for zooming in [x1 x2 y1 y2] 
%             [default is whole frame]
%    refidx = frame index. If not input or empty [def]: calc difference from current frame
%             If scalar input: calc difference from this single reference frame
%
% OPTIONAL OUT:  Fmean = average pixel difference between frames (Abs value)
%                Fmax  = maximum pixel difference between frames (Abs value)
%                idx   = vector of frame index - centered on the half value
%                        between differenced frames (e.g. [1.5 2.5 ... N-0.5]

% C Rowell March 2020
%
%
clear textprogressbar
nargoutchk(0,3)

    if nargin<6
        refidx = [];
    end
    if nargin<5
        roi = [];
    end
    if nargin<4
        regHeads = [];
    end
    if nargin<3
        regDir = [];
    end
    
    if ~isempty(regDir)
        compareReg = true;
        if isempty(regHeads)
            regHeads = fullfile(regDir,'frameHeads.mat');
        end
        if ~exist(regHeads,'file')
            error('Cannot find registered parameter table')
        else
            load(regHeads)
            TR = T;
        end
    else
        compareReg = false;
    end

    load(matHeads)

    load(fullfile(matDir,T.File{1}))
    imsz = size(Frame);
    if isempty(roi)
        roi1 = [1 imsz(2) 1 imsz(1)];
        load(fullfile(regDir,TR.File{1}))
        imsz2 = size(Frame);
        roi2 = [1 imsz2(2) 1 imsz2(1)];
    else
        roi1 = roi;
        roi2 = roi;
    end
    
    
    if isempty(refidx)
        use_ref = false;
    else
        use_ref = true;
    end

    if compareReg
        idx = {cellfun(@str2num, T.Properties.RowNames)', cellfun(@str2num, TR.Properties.RowNames)'};
        Fmean = {zeros(length(idx{1})-1,1) zeros(length(idx{2})-1,1)};
        didx{1} = idx{1}(1:end-1)+0.5;
        didx{2} = idx{2}(1:end-1)+0.5;
    else
        idx = {cellfun(@str2num, T.Properties.RowNames)'};
        Fmean = {zeros(length(idx{1})-1,1)};
        didx{1} = idx{1}(1:end-1)+0.5;
    end
    Fmax  = Fmean;
%     didx = idx(1:end-1)+0.5;
    
%     figure
    textprogressbar('MAT frame difference plots > ')
    for ix = 1:length(idx{1})-1
        
        
        if ix==1
            if use_ref
                load(fullfile(matDir,T.File{num2str(refidx)}))
                F1 = Frame(roi1(3):roi1(4),roi1(1):roi1(2));
            else
                load(fullfile(matDir,T.File{num2str(idx{1}(ix))}))
                F1 = Frame(roi1(3):roi1(4),roi1(1):roi1(2));
            end
        elseif ~use_ref
            F1 = F2;
        end
        load(fullfile(matDir,T.File{num2str(idx{1}(ix+1))}))
        F2 = Frame(roi1(3):roi1(4),roi1(1):roi1(2));
        
        Fmean{1}(ix) = mean(abs(F2(:)-F1(:)));
        Fmax{1}(ix)  = max(abs(F2(:)-F1(:)));
        
        textprogressbar(ix/length(idx{1})*100)
    end
    textprogressbar(' -> Done')

    if compareReg
%         idx{2} = cellfun(@str2num, TR.Properties.RowNames)';
%         Fmeanr = zeros(length(idx{1})-1,1);
%         Fmaxr  = Fmeanr;
%         didx{2} = idx{2}(1:end-1)+0.5;

    %     figure
        textprogressbar('Registered frame difference plots > ')
        for ix = 1:length(idx{2})-1


            if ix==1
                if use_ref
                    load(fullfile(regDir,TR.File{num2str(refidx)}))
                    F1 = Frame(roi2(3):roi2(4),roi2(1):roi2(2));
                else
                    load(fullfile(regDir,TR.File{num2str(idx{2}(ix))}))
                    F1 = Frame(roi2(3):roi2(4),roi2(1):roi2(2));
                end
            elseif ~use_ref
                F1 = F2;
            end
            load(fullfile(regDir,TR.File{num2str(idx{2}(ix+1))}))
            F2 = Frame(roi2(3):roi2(4),roi2(1):roi2(2));

            Fmean{2}(ix) = mean(abs(F2(:)-F1(:)));
            Fmax{2}(ix)  = max(abs(F2(:)-F1(:)));

            textprogressbar(ix/length(idx{1})*100)
        end
        textprogressbar(' -> Done')

    end
    
    figure
    ax(1)=tightSubplot(2,1,1,[],0);
    plot(didx{1},Fmean{1},'LineWidth',1.2)
    if compareReg
        hold on
        plot(didx{2},Fmean{2},'LineWidth',1.2)
        legend('Raw frames','Registered frames')
    end
    if use_ref
        ylabel(sprintf('Mean(|T_{i+1} - T_{%i}|) [K]',refidx))
    else
        ylabel('Mean(|T_{i+1} - T_i|) [K]')
    end
    
    ax(2)=tightSubplot(2,1,2,[],0);
    plot(didx{1},Fmax{1},'LineWidth',1.2)
    if compareReg
        hold on
        plot(didx{2},Fmax{2},'LineWidth',1.2)
    end
    xlabel('Frame Index Pair')
    if use_ref
        ylabel(sprintf('Max(|F_{i+1} - F_{%i}|) [K]',refidx))
    else
        ylabel('Max(|T_{i+1} - T_i|) [K]')
    end
    linkaxes(ax,'x')
%     if use_ref
%         legend(sprintf('Mean(|F_{i+1} - F_{%i}|)',refidx),sprintf('Max(|F_{i+1} - F_{%i}|)',refidx))
%     else
%         legend('Mean(|F_{i+1} - F_i|)','Max(|F_{i+1} - F_i|)')
%     end
if nargout>=1
    varargout{1} = Fmean;
end
if nargout>=2
    varargout{2} = Fmax;
end
if nargout==3
    varargout{3} = didx;
end

end