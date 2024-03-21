function varargout = irbAsc2Mat(idir,odir,param_file,glob_spec)
% irbAsc2Mat(idir,odir,param_file,glob_spec)
% Load IRBIS ascii and convert to .mat files compatible with plumeTracker
%
% IN:   idir = directory containing input ascii images. File names should
%               be formatted numerically and separated by an underscore. 
%               e.g. file_001.txt, file_002.txt, etc
%       odir = directory to save output mat files
%
% OPTIONAL IN:
%       param_file = text file path containing meta data from IRBIS files
%                    Defaults to empty [], uses individual file metadata
%                    instead
%                   ALTERNATE: enter as 'None' to skip writing a param
%                   file (for use with irbHeads2Table)
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
%  C Rowell, August 2018
%
% ------- irbAsc2Mat example input: -------
% 
% datadir   = '/myhome/mydata/';
% ascDir  = fullfile(datadir,'ascii/');
% matDir  = fullfile(datadir,'mat/');
% 
% wdir = fullfile(homedir,'Kahuna/data/sabancaya_5_2018/image_exports/test_data/');
% idir = ascDir;
% odir = matDir;
% 
% asc_params = fullfile(datadir,'BI052500_conv_1807.txt');
%   OR
% asc_params = 'None'; % if converted file from 'irbHeads2Table.m'
% 
% glob_spec = '.txt';
%
% [T,IRTheads_file] = irbAsc2Mat(ascDir,matDir,asc_params',glob_spec);

fprintf('\n========= IRBIS Ascii to Mat Conversion =========\n')

% ============================= DO THE THING =============================
if nargin<3
    param_file = [];
end
if nargin<4
    glob_spec = '*.txt';
end
if isempty(glob_spec)
    glob_spec = '*.txt';
end
if strcmp(param_file,'None')
    writeParams = false;
else
    writeParams = true;
end
nargoutchk(0,2);

disp('Converting ASCII to MAT files...')
fprintf('  In  dir:\t%s\n',idir)
fprintf('  Out dir:\t%s\n',odir)
fprintf('  Glob specifier:\t%s\n',glob_spec)
if ~exist(odir,'dir')
    fprintf('Making new output directory.')
    mkdir(odir)
end

% Get file list from directory, S, and exported parameter table, T
Slist = glob(fullfile(idir,glob_spec));

if ~isempty(param_file)
    % Cut out param file if it's in the glob
    [aa,ii] = ismember(param_file,Slist);
    if aa
        Slist(ii)=[];
    end
end

[paths,names,exts] = cellfun(@fileparts,Slist,'UniformOutput',false);
% digs = numel(num2str(size(S,1))); % number of digits for file indices
% Sidx = cellfun(@(x) str2num(x(end-digs+1:end)) ,names);
if isempty(Slist)
    error('Ascii file list is empty! Check your input directory!')
end
noms = splitfix(names,'_');
Sidx = double(noms(:,end));


% Load param file if it exists
if and(~isempty(param_file),writeParams)
    T = readIRBISparams(param_file);
    % Use glob spec to cut down param table
    Tnew = T(cellstr(string(Sidx)),:);
    use_table = true;
    
    fprintf('Parameter file lists %i images, glob found %i files\n',size(T,1),numel(Slist))
else
    fprintf('No parameter file, Glob found %i files\n',numel(Slist))
    use_table = false;
end

% Loop through file list
% fprintf('Converting %i/%i files...\n',numel(S),size(T,1))

% GOOD SO FAR TO HERE
disp('Reading ascii files...')
for ss = 1:length(Slist)
    % Get image index
    name = names{ss};
%     [path,name,fext] = fileparts(S{ss});
    idx = Sidx(ss); %double(string(name(end-2:end)));
    
    % Assuming data in digital format for now
%% Table management
    if use_table
%         Frame = double(imread(fullfile(idir,T.File{idx})));
%         digs = numel(num2str(size(T,1)));
        Frame = readIrbAsciiFrame(Slist{ss});
%         Frame = readIrbAsciiFrame(strrep(T.File{idx},'.irb','.txt'));
        T.File{idx} = [name '.mat'];
        oidx = idx;
        
    else
        Frame = readIrbAsciiFrame(Slist{ss});
   
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

    %             lst2 = struct2cell(lst);
    %             f = fieldnames(lst);
    %             lst2 = [lst2; num2cell(nan(numel(a(a)),1)) ];
    %             f = [f; T.Properties.VariableNames(a)'];
    %             lst = cell2struct(lst2,f);
            end
            T(ss,:) = struct2table(lst);
        end
        oidx = ss; % Bad workaround, but hey
    end
    

    File_DateTime = T.Timestamp(oidx);
    Oidx(ss) = idx;

    %% Save and output
    fprintf('Saving: %i,\t%s\n',idx,T.File{oidx})
    save(fullfile(odir,T.File{oidx}),'Frame','File_DateTime')
end

% Fix table output as in readIRBISParams if no param file
if ~use_table
    T.Properties.RowNames = cellstr(string(Oidx));
    T = fixParamTable(T);
else
    T = T(Oidx,:); % A bit unneccesary without a param file
end

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

ofile = fullfile(odir,'headsIRBASC2MAT.mat');

% Pull param file data only for grabbed frames
if writeParams
    save(ofile,'T')
else
    disp('NO parameter file written.')
end

if nargout>=1
    varargout{1} = T;
end
if nargout==2
    varargout{2} = ofile;
end
disp('Done!')
end

function [heads] = readIrbAsciiHead(fname)
%    Read it in, mofugga
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

function [frame] = readIrbAsciiFrame(ifile)
%    Read it in, sucka
%   OUT:    Frame = image array

% Reference struct for corresponding header names

    % Read image
    frame = readtable(ifile,'delimiter','tab','ReadVariableNames',false);
    frame = double(string(strrep(table2cell(frame),',','.')));
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