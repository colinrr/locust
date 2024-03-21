function varargout = irbHeads2Table(idir,odir,glob_spec)
% irbHeads2Params(idir,odir,glob_spec)
% Load IRBIS ascii, read headers only, and convert to .mat parameters file
% compatible with plumeTracker
%
% IN:   idir = directory containing input ascii images. File names should
%               be formatted numerically and separated by an underscore. 
%               e.g. file_001.txt, file_002.txt, etc
%       odir = directory to save output mat files
%   OPTIONAL IN:
%       glob_spec  = GLOB expression to grab a subset of files 
%                   (eg. '*_?23.tif') [Default is '*.txt']
%
% OUT: 
%  ->  .mat files saved to directory "odir", along with header file
%      'headsIRBASC2MAT.mat' IF 'param_file' was specified
%
%   varargout{1} = T : table of frame header information
%   varargout{2} = ofile : path to saved output header table
%
%  C Rowell, March 2020

% ------- irbHeads2Table example input: -------
% 
% datadir   = '/myhome/mydata/';
% ascDir  = fullfile(datadir,'ascii/');
% matDir  = fullfile(datadir,'mat/');
% 
% glob_spec = '.txt';
%
% [T,heads_file] = irbHeads2Table(ascDir,matDir,glob_spec);

    fprintf('\n========= Convert IRBIS Ascii Heads to Param File =========\n')

    % ============================= DO THE THING =============================
    if nargin<3
        glob_spec = '*.txt';
    end
    if isempty(glob_spec)
        glob_spec = '*.txt';
    end

    fprintf('  In  dir:\t%s\n',idir)
    fprintf('  Out dir:\t%s\n',odir)
    fprintf('  Glob specifier:\t%s\n',glob_spec)
    if ~exist(odir,'dir')
        fprintf('Making new output directory.')
        mkdir(odir)
    end

    % Get file list from directory, S, and exported parameter table, T
    Slist = glob(fullfile(idir,glob_spec));

    [paths,names,exts] = cellfun(@fileparts,Slist,'UniformOutput',false);
    % digs = numel(num2str(size(S,1))); % number of digits for file indices
    % Sidx = cellfun(@(x) str2num(x(end-digs+1:end)) ,names);
    if isempty(Slist)
        error('Ascii file list is empty! Check your input directory!')
    else
        fprintf('  Generating parameter table, Glob found %i files\n',numel(Slist))
    end
    noms = splitfix(names,'_');
    Sidx = double(noms(:,end));

    %%
    textprogressbar('  Reading headers: ')
    for ss = 1:length(Slist)
        % Get image index
        name = names{ss};
        idx = Sidx(ss); 

        head  = readIrbAsciiHead(Slist{ss});

        % Initialize table
        myvars = genvarname(strrep(head(:,1),'-','_'));
        lst    = cell2struct(head(:,2),myvars);
        lst.File = {[name '.mat']};
        if ss==1
    %             T = cell2table(head(:,2)','VariableNames',myvars);
            T = struct2table(lst);
        else

            % Add entries for missing values
            a = ~ismember(T.Properties.VariableNames,fieldnames(lst));
            if any(a)
                b = find(a);
                for jj = 1:length(b)
                    ii = b(jj);
                    varname = T.Properties.VariableNames{ii};
                    varval  = T{end,ii};
                    if isnumeric(varval)
                        lst.(varname) = NaN;
                    elseif iscell(varval)
                        lst.(varname) = {''};
                    end
                end
            end
            T(ss,:) = struct2table(lst);
        end
        oidx = ss; % Bad workaround, but hey

        File_DateTime = T.Timestamp(oidx);
        Oidx(ss) = idx;
        textprogressbar(ss/length(Slist)*100)
    end
    % Fix table output as in readIRBISParams if no param file
    T.Properties.RowNames = cellstr(string(Oidx));
    T = fixParamTable(T);

    % Some checks
    [~,Ti1] = sort(str2double(T.Properties.RowNames));
    T = T(Ti1,:);
    % [T,Ti1] = sortrows(T,'RowNames');
    [~,Ti2] = sortrows(T,'Time');
    if ~isequal(Ti1,Ti2)
        fprintf('\n\n')
        warning('Sorting frame table by Index does not match sorting table by Time. TIME VECTOR OR INDICES MAY BE INCORRECT!')
        fprintf('\n')
    end
    if numel(T.Time)~=numel(unique(T.Time))
        fprintf('\n\n')
        warning('Repeat values found in time vector! TIME VECTOR OR INDICES MAY BE INCORRECT!')
        fprintf('\n')
    end

    oname = fullfile(odir,'frameHeads.mat');
    if nargout>=1
        varargout{1} = T;
    end
    if nargout==2
        varargout{2} = oname;
    end
    % Pull param file data only for grabbed frames
    save(oname,'T')
    textprogressbar(' -> Done!')
end

function [heads] = readIrbAsciiHead(fname)
%    Read it in
%   OUT:    head = head information in row vector, ready for table.

% Reference struct for corresponding header names
hnames = struct('Name','File');
heads = {};
    % Read header
    fid = fopen(fname);
    tline = fgetl(fid); % Skip writing the first line
    while ~strcmp(tline,'')
        tline = fgetl(fid);
        chunks = splitfix(tline,':');
        varname = chunks(1);
        
        % Parse header lines a bit
        if ~strcmp(tline,'')
            if strcmp(chunks(1),'Name')
                varname = 'File';
                ff = splitfix(chunks(end),{'\','/'});
                varval = {char(ff(end))};
            elseif strcmp(chunks(1),'Date')
                varname = 'Timestamp';

                vv = splitfix(tline,{'.',':',' '});
                if iscell(vv); vv = string(vv); end
                D = fliplr(vv(2:4)'); T = vv(5:end)';
                varval = datenum(double([D T]));
            elseif strcmp(chunks(1),'ms')
                varname = 'msec';
                varval = double(strrep(chunks(end),',','.'));
            else
                varname = char(chunks(1));
                varval = double(strrep(chunks(end),',','.'));
                if isnan(varval)
%                     varval = char(join(chunks(2:end)));
                    varval = cellstr(string(join(chunks(2:end))));
                end
            end
            heads = [heads; {varname} {varval}];
        end
    end
    fclose(fid);

end

function T = fixParamTable(T)
% For cases where no IRBIS param file is present, fix the table output (eg
% timestamps and indices) in the table derived from individual files

    % Use timestamps in case no msec information (and because bullshit IRB
    % files don't record absolute start time to greater than 1 s precision)
    Ms   = T.msec;
    timestamp = T.Timestamp; %datenum([fliplr(Date) Time(:,1:2) Time(:,3)]);
    t0 = timestamp(1);
    reltime = Ms/1000;
    reltime = reltime-reltime(1); % Zero out reltime in case not starting at beginning of irb file

    rel_timestamp = timestamp-t0;
    
    % ----------------------------------------------------------------
    %  THIS BIT fixes (attempts to) issues where the image sequence
    %  contains both time-lapse and real-time frames, as the time-lapse
    %  frames in .IRB files take time stamps only to the nearest second
    %  (stupid, I know) and so cannot be exactly aligned with the
    %  millisecond values in real-time images. The resulting timing error
    % is of order 1 second
    if ~any(reltime>0) % Quick fix for only time lapse imagery
        reltime = rel_timestamp*86400;
    elseif sum(reltime==0)>1 % Fix for mixed frames
        % Bullshit approach to fix added frames at start
        bsIdx = find(reltime~=0,1,'first'); % Index of first
        fudge_vec = ones(size(rel_timestamp))*rel_timestamp(bsIdx-1)*86400;
        fudge_vec(1:bsIdx-1) = rel_timestamp(1:bsIdx-1)*86400;
        reltime = reltime+fudge_vec;
    end
    % Otherwise (only real-time frames), the reltime vector should be fine
    % ------------------------------------------------------------------
    
    % Bullshit approach to fix added frames at start
%     bsIdx = find(reltime~=0,1,'first'); % Index of first
%     fudge_vec = ones(size(rel_timestamp))*rel_timestamp(bsIdx-1)*86400;
%     fudge_vec(1:bsIdx-1) = rel_timestamp(1:bsIdx-1)*86400;
%     reltime = reltime+fudge_vec;

    T.Time = reltime;
    T = [T(:,1:2) T(:,end) T(:,3:end-1)]; % Rearrange reltime to front of table
    
end

function st = splitfix(st,stformat)
% Attempsts to fix a weird issue wheres split outputs either cells or
% strings depending on the machine/matlab version
    st = split(st,stformat);
    if iscell(st)
        st = string(st);
    end
end