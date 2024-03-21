function [T,idxCut] = fixTableTimes(paramIRB, paramIRT,matDir)
% Function to load frame tables (output from irbAsc2Param.m AND/OR
% irbHeads2Table ) and fix time vectors by finding dropped frames and
% loading accordingly.
%
% ISSUE: Correct (to ms) time stamps come from straight conversion from raw
% .irb videos, but these do not output correct data. So we need a separate
% conversion of IRT_analyzer, which only outputs timestamps to the second,
% BUT also drops some frames for no obvious reason. So the challenge is to
% figure out WHICH frames were dropped, without the benefit of timestamps
% with greater than 1 second accuracy.
%
% IN:   paramIRB = param file (Table) from .IRB conversion. Correct
%               timestamps, ALL frames
%       paramIRT = param file (Table) from .IRT conversion. Timestamps to
%               nearest second only, some dropped frames
%       matDir   = directory of .mat frames. Function uses image
%               differencing to find dropped frames
%
% OUT: T = updated table with complete millisecond time, absent dropped frames
%      idxCut = index location of dropped frames

    fprintf('\n========= Find dropped frames, fix header tables =========\n')
    fprintf('Reference headers:\n\t%s\n',paramIRB)
    fprintf('Data headers:\n\t%s\n',paramIRB)

    load(paramIRB)
    Tirb = T;
    load(paramIRT)
    Tirt = T;

    Nirb = size(Tirb,1);
    Nirt = size(Tirt,1);
    Ndrop = Nirb-Nirt;
    
    fprintf('Total frames:\t%i\nIRT frames:\t%i\nDropped frames:\t%i\n',Nirb,Nirt,Ndrop)

    ts0 = min(Tirb.Timestamp);
    ms = Tirt.msec;
    ts = (Tirt.Timestamp-ts0)*86400;
    idx = str2double(Tirt.Properties.RowNames);

    ms_raw = Tirb.msec;
    ts_raw = (Tirb.Timestamp-ts0)*86400;
    idx_raw = str2double(Tirb.Properties.RowNames);
    ms_out = ms_raw; % time series for editing into IRT table
    idx_edit = idx_raw; % Not output, but useful for checking
    ts_edit = ts_raw;

    
    % Get first switch between seconds
    switchI = find(Tirb.Timestamp~=ts0,1,'first');
    % USE IRT output indices
    % Loop over these, find when ts does not match, adjust ts_raw to fit by
    % killing that index,ts, and ms value. Then add correct ms value to IRT
    % table and continue.
    %%
    idxCut = zeros(Ndrop,1);
    dropCount = 0;
    for ii = 1:Nirt
        % Trigger when timestamps become unequal or check last samples
        if ts(ii)~=ts_edit(ii) ||  and(ii==Nirt,length(ts_edit~=Nirt))
            count = 0;
            trigt = ts(ii);
            % Find number of consecutive frames that are unequal
            while trigt ~= ts_edit(ii+count)
                count = count+1;
                trigt = ts(ii+count);
            end
            % Get indices of last second
            idxCheck = find(ts==ts(ii-1));
            idxComp   = find(ts_edit==ts(ii-1)); % Indices to compare against
            
            % Check for final second
            if and(ii==Nirt,length(ts_edit~=Nirt)) && count==0
                count = length(idxComp)-length(idxCheck);
            end
            % Load frames and find differences
            dF = zeros(length(idxCheck)-1,1);
            for ix = 1:length(idxCheck)-1
                if ix==1
                    load(fullfile(matDir,Tirt.File{num2str(idxCheck(ix))}))
                    F1 = Frame;
                else
                    F1 = F2;
                end
                load(fullfile(matDir,Tirt.File{num2str(idxCheck(ix+1))}))
                F2 = Frame;
                % Simple check for larger changes
                dF(ix) = max(abs(F2(:)-F1(:)));
            end
            % Get <count> largest dF values
            [dFsort,Isort] = sort(dF,'descend');
            Isort = Isort(1:count) + 1; %-> dF is large when NEXT frame is missing
            % Clear ms value(s), update time series for iterating
            dropCount = dropCount+length(Isort);
            idxCut(dropCount:dropCount+length(Isort)-1) = idxComp(Isort);
            ms_out(idxComp(Isort)) = []; % Clip missing timestamps
            ts_edit(idxComp(Isort)) = [];
            idx_edit(end-(count-1):end) = []; % Cut indices from end
        end
    end
    
    % CHECK
    if dropCount>Ndrop
        error('Too many frames dropped...')
    elseif dropCount<Ndrop
        if ~any(ts_edit(1:length(ts))~=ts)
            warning('Dropping missing frames from end of sequence.')
            idxCut(dropCount+1:end) = length(ts)+1:length(ts_edit);
            dropCount = dropCount+(length(ms_out)-length(ts));
            ms_out = ms_out(1:length(ts));
        else
            error('Could not find all dropped frames - check time series manually!')
        end
    end
    
    % Update correct timestamps
    Tirt.Time = ms_out/1000;
    Tirt.msec = ms_out;
    
    % Do a check
    [~,Ti1] = sort(str2double(Tirt.Properties.RowNames));
    T = Tirt(Ti1,:);
    % [T,Ti1] = sortrows(T,'RowNames');
    [~,Ti2] = sortrows(T,'Time');
    check1 = ~isequal(Ti1,Ti2);
    check2 = numel(T.Time)~=numel(unique(T.Time));
    if check1
        fprintf('\n\n')
        warning('Sorting frame table by Index does not match sorting table by Time. TIME VECTOR OR INDICES MAY BE INCORRECT!')
        fprintf('\n')
    end
    if check2
        fprintf('\n\n')
        warning('Repeat values found in time vector! TIME VECTOR OR INDICES MAY BE INCORRECT!')
        fprintf('\n')
    end
    if and(~check1,~check2)
        ofile = fullfile(matDir,'frameHeadsFixed.mat');
        fprintf('Writing new parameter file:\n\t%s\n',ofile)
        save(fullfile(matDir,'frameHeadsFixed.mat'),'T')
    else
        disp('Still some timing issues, so skipped writing a new parameter file.')
    end
end