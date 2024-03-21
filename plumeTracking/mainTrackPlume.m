function [Tout,Ref ] = mainTrackPlume(inputDir,heads,outputDir,ref,deb,fin,dN,fixZ,plot_flags,cax )
% function [Tout,Ref ] = mainTrackPlume(inputDir,outputDir,ref,deb,fin,dN,fixZ,plot_flags )
% MAINTRACKPLUME Main routine for plumeTracker algorithm (Bombrun et al., 2018).
% REQUIRED INPUT:    
%   inputDir = directory of .mat files w/ image data
%               - optionally should contain file params.mat with IRBIS frame meta data
%               - files should be indexed with integer
%                 names, eg file_001.mat, file_002.mat, etc
%               - should contain variables:
%                    > Frame: image array
%                    > File_DateTime: timestamp of image in
%                            vector form: [Y Mo D H Mi S ms]
%                   
% OPTIONAL INPUT:\
%   heads      = header table. [default search: inputDir/frameHeads.mat]
%   outputDir  = [default: inputDir/results/ ]
%   ref        = reference image for background removal. Normally the same as deb
%   deb        = starting image [default: minimum file index found]
%                 !! ALTERNATIVE: deb may be specified as strictly
%                 increasing vector of indices, which will be used to
%                 compile the image list INSTEAD of fin and dN. Ref will
%                 remain unchanged.!!
%   fin        = final image [default: max file index found]
%   dN         = image step - may need to jump images to get sufficient
%                differences between frames - masks can be interpolated to
%                skipped frames later
%   fixZ       = 2 element vector [x y] (pix) giving a fixed position for the plume base
%                (prevents segmentation of anything below this position
%   plot_flags = boolean vector to save image output: [plume_params png gif video]
%                   -> plume_params tell whether to plot measured plume
%                   parameters on the images
%                   -> default is [false false false true] (aka output
%                   video only)
%   cax        = temperature range for plotting
%
% OUTPUT:   
%   Tout = updated image table with masks, positions, etc
%   Ref  = structure of plumeTracking parameters

% Original segmentation algorithm written by Maxime Bombrun
% -> Bombrun, M., Jessop, D., Harris, A., & Barra, V. (2018). 
% -> An algorithm for the detection and characterisation of volcanic plumes
% -> using thermal camera imagery. Journal of Volcanology and Geothermal 
% -> Research, 352, 26?37. https://doi.org/10.1016/j.jvolgeores.2018.01.006

% mainTrackPlume.m and trackingPlume.m modified by Colin Rowell, July 2018,
% to streamline I/O for use with IRBIS data

% C.Rowell CHANGES:
% DONE  1) generalized file indexing to use image header tables
% DONE  2) change image loading to target the specific file name 
% DONE  3) Incorporate image timestamp information
% DONE  4) New function I/O wrapper
% DONE  5) Update data table output
% DONE    % 5a) Set data to use input data table from IRBIS
% DONE    % 5b) Output similar table
%         % 5c) Show ref image in small screen
% DONE  6) Incorporate output mask and contour data into tables
% DONE  7) Output extracted mask widths, and heights, and positions
% DONE  8) Optionally fix plume base to a specific pixel (fixZ)
% DONE  9) Set plot flags (eg gif, video) as optional input arguments

fprintf('\n========= plumeTracker =========\n')

%% Sample input: Load data set and Initialisation
% inputDir='/home/crowell/Kahuna/data/plumeTracking/data/stromb/';
% outputDir='/home/crowell/Kahuna/data/plumeTracking/data/stromb/results/';
% inputName  = '/home/crowell/Kahuna/data/sabancaya_5_2018/image_exports/BI0525_big_explosion/raw_values/mat/';
% outputName = fullfile(inputName,'PTresults/');

% ref = 93; % Index of reference image to use
% deb = 93; %? Need not be the same as reference image
% fin = []; % Leave empty to grab all files to end
% freq=15;
% dt = 5;

% Plotting
% write_video = true;
% write_gif   = true;
%% Input parsing
narginchk(1,10)
if nargin<2; heads = []; end
if nargin<3; outputDir = fullfile(inputDir,'PTresults/'); end
if nargin<4; ref = []; end
if nargin<5; deb = []; end
if nargin<6; fin = []; end
if nargin<7; dN  = 1;  end
if nargin<8; fixZ = []; end
if nargin<9; plot_flags = [false false false true]; end
if nargin<10; cax=[]; end

plume_params = plot_flags(1);
png = plot_flags(2);
gif = plot_flags(3);
vid = plot_flags(4);

if ~exist(inputDir,'dir')
    error(sprintf('Input directory does not exist:\n\t%s',inputDir))
end

fprintf('Plume Tracking in directory:\n\t%s\n',inputDir)
if isempty(heads)
    heads = fullfile(inputDir,'frameHeads.mat');
end
if isempty(outputDir)
    outputDir = fullfile(inputDir,'results/');
end
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

% Load file list
[listImg,nbFrame,refS,deb,fin]=parseInputDir(inputDir,heads,deb,ref,fin,dN);
ref = refS.file;
if ~any(listImg.Time>0) % Quick fix for timelapse files
    tvec = (listImg.Timestamp - min(listImg.Timestamp))*86400;
else
    tvec = listImg.Time - min(listImg.Time);
end
[~,ilen,~] = fileparts(listImg.File{1});
% Get number of index digits
ilen = split(ilen,'_');
ilen = ilen(end);
ilen = numel(ilen{1});
%% Load reference image and compute time step
fprintf('Reference Image: %s\n',ref)
[src,tr]=loadImg(fullfile(inputDir,ref)); 

% Computing time step just once has problems for any data set where the
% level of activity changes or time steps are not consistent...
% dt=computeTimeStep(inputName, nbFrame);

noutFrame = nbFrame-1; %length(1+dN:dN:nbFrame);
% Frame - index, name, timestamp/rel time, height, width, h/w positions,
% angles?
entete={'Index' 'Frame' 'VidTime' 'PixHeight' 'PixWidth' 'Outline' 'Mask' 'Positions'};
ent_units = {'' '' '(s)' '(pix)' '(pix)' 'bool' 'bool' '(pix)'};
ent_desc  = {'File index', 'Frame number in video, etc','Relative time from start frame',...
    'Plume height in pixels','Plume width in pixels',...
    'Outline from plumeTrack', 'Mask from plumeTrack',...
    'Position of plume scale measurements: Rows= [pt1, pt2],  Columns= [WidthX, HeightX, WidthZ, HeightZ]'};
contenu=cell(noutFrame,8);

[~,zInit]=max(src(:));
[zInit,~]=ind2sub(size(src),zInit);

maskTSave=zeros(size(src,1),size(src,2),10);
idx=1;
%% %%%%%%%%%%%%%%%% VIDEO and PLOT SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open a video writer
if vid
    vidObj = VideoWriter(fullfile(outputDir, 'mapT'),'Motion JPEG AVI');
    vidObj.FrameRate = 10;
    open(vidObj);
end
% Open a figure and parameterize it
fig=figure;
set(fig, 'Position', [100 100 1*size(src,2) 1*size(src,1)])
axis([0 size(src,1) 0 size(src,2)]);
set(gcf, 'PaperPositionMode', 'auto');
set(gca,'position',[0 0 1 1],'units','normalized','XColor',[0.9 0.9 0.9])
% colormap(jet);

%% ================ Run through all images ================
fprintf('Start/stop indices: deb = %i, fin = %i\n',deb,fin)
% fprintf('Found %i file indices in range: %s to %s\n',nbFrame, listImg.Properties.RowNames{1},listImg.Properties.RowNames{end})
fprintf('Found and Tracking %i images, with indices %s to %s (used step dN = %i)...\n', nbFrame, listImg.Properties.RowNames{1},listImg.Properties.RowNames{end},dN)
for i=1:nbFrame-1 %dN:nbFrame-dN
% for i=deb: dt:fin-dt
    % Extract mask

%     fprintf('%i: %s\n',listImg{i+dt,1},inputDir)
    Fidx = str2num(listImg.Properties.RowNames{i+1});
    [currSrc,ts]=loadImg(fullfile(inputDir,listImg.File{i+1}));
    if i==1; t0 = ts; end
    fprintf('Frame %i) file %i:\t%s\n',i,Fidx,datestr(ts,'yyyy-mm-dd HH:MM:SS.FFF'))
    [outline, mask] = trackingPlume(i,ref,inputDir,listImg,dN);
    
    
    zoneT=find(mask==1);
    [z0,x0]=ind2sub(size(mask), zoneT);
    
%     abs_t = datenum(ts);
%     rel_t = datenum(ts)-datenum(t0);
    rel_t = tvec(i+1);
   %% Extract plume parameters
   % - Heigth, Width, ?Velocity?
    [Pt,pti]=min(z0); % Plume top?
%     if zInit ~= max(z0)
    if isempty(fixZ)
        [zInit,zi]=max(z0);
        xInit = x0(zi);
    else
        xInit = fixZ(1);
        zInit = fixZ(2);
    end
%     end
    plumeHeight=abs(zInit-Pt);          % Absolute max 
%     plumeWidth=abs(max(x0)-min(x0));    % Absolute maximum width
    [plumeWidth,zW]=max(sum(mask,2));   % Max width in any one horizontal line - crude atm
    xW0 = find(mask(zW,:),1,'first');
    xW1 = find(mask(zW,:),1,'last');
    % Rows: pt1, pt2.  Columns: WidthX, HeightX, WidthZ, HeightZ
%     hwind = [x0(pti) x0(zi) zInit Pt ; xW0 xW1 zW zW]; 
%     hwind = [xW0 x0(pti) zW Pt; xW1 xInit zW zInit];
    %        wx1 wy1 wx2 wy2 hx1 hy1 hx2 hy2]
    hwind = [xW0 zW xW1 zW x0(pti) Pt xInit zInit];
    %% Save the images
    imagesc(mask)
%     saveas(gcf,['./' outputName '/mask_' num2str(i+step,'%0.4d') '.png'],'png');
    if png
        saveas(gcf,fullfile(outputDir, ['mask_' num2str(Fidx,['%0.' num2str(ilen) 'd']) '.png']),'png');
    end
    plotContourImage(currSrc,outline,cax) % quick imagesc setup
%     imagesc(outline)
%     caxis(cax)
    colorbar('east','Color','w')
     
    if png
%     saveas(gcf,['./' outputName '/contour_' num2str(i+step,'%0.4d') '.png'],'png');
    saveas(gcf,fullfile(outputDir, ['contour_' num2str(Fidx,['%0.' num2str(ilen) 'd']) '.png']),'png');
    end
    
%     mapT=currSrc.*mask;
    % fin=deb+nbFrame-1;

    %% Set the Gif
    maskT=mask;
    for z=Pt:zInit
        slice=mask(z,:).*currSrc(z,:);
        list=slice~=0;
        M=mean(slice(list));
        maskT(z,:)=M.*list;
    end
    maskTSave(:,:,idx)=maskT;
    
    %%%%%%%%%%%%%%%%%% VIDEO %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    plotContourImage(currSrc,outline,cax)
%     imagesc(outline);
%     caxis(cax)
    %     caxis([250 cax(2)]) % TEMPORARY
    colorbar('east','Color','w')

    
    if plume_params
        plotHW(hwind)
    end
%     timing= (i+dN)/30;
%     M=floor(timing/60);
%     S=num2str(floor(timing-M*60));
%     M=num2str(M);
%     timing= [sprintf('%02s',M) ':' sprintf('%02s',S)];
    leg=sprintf('%i\n%s',Fidx,datestr(rel_t/86400,'HH:MM:SS.FFF'));
    t=text(0.1,0.85,leg,'FontSize',14,'Color', [1 1 1],'Units','normalized');
    %scale if known
%     rectangle('Position', [size(src,2)-50, size(src,1)-50, 4, 34],'FaceColor', [1 1 1])
%     t2=text(260,217,'200 m','FontSize',20,'Color', [1 1 1]);
    %colorbar
    hold off;
    axis off
%     set(gca,'XColor',[0.9 0.9 0.9])
    % Transform the frame in video
        F = getframe(fig);
    if vid
        writeVideo(vidObj,F);
    end
    %%%%%%%%%%%%%%%%%% VIDEO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% GIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if gif
        im = frame2im(F);
        [imind,cm] = rgb2ind(im,256);
        if i == 1
            imwrite(imind,cm,[outputDir 'mapT.gif'],'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,[outputDir 'mapT.gif'],'gif','WriteMode','append');
        end
    end
    %%%%%%%%%%%%%%%%%%% GIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Save information in 'content'
%     time=i/freq;
%     [H,M,S]=secs2hms(time);
%     time=[sprintf('%02.0f', H) ':' sprintf('%02.0f', M) ':' sprintf('%02.0f', S)];
%     time = rel_t; %datestr(rel_t, 'HH:MM:SS.FFF');
    if i==1
            %        wx1 wy1 wx2 wy2 hx1 hy1 hx2 hy2]
        hwind = [xW0 zW xW1 zW x0(pti) Pt xInit zInit];

        contenu(1,:)={listImg.Properties.RowNames{i} 0 tvec(i) 0 0 [] sparse(mask.*0) repmat(fixZ,[1 4])};
    end
    contenu(idx+1,:)={listImg.Properties.RowNames{i+1} idx rel_t plumeHeight plumeWidth outline sparse(mask) hwind};
    idx=idx+1;
end

% Close the objects
% if exist('fig')
%     close(fig);
% end
if exist('vidObj','var')
    close(vidObj);
end

%% Write table and Save the content
contenu = cell2table(contenu,'VariableNames',entete);
contenu.Properties.VariableDescriptions = ent_desc;
contenu.Properties.VariableUnits = ent_units;
tt = listImg(contenu.Index,:);
Tout = [tt(:,1) contenu(:,2) tt(:,2) contenu(:,3:end) tt(:,3:end)];

Ref = struct('idir',inputDir,'odir',outputDir,'refFile',ref,'refIm',src,'refParams',refS.par);

output_time = datestr(now);
if ~isempty(contenu)
%     filename=fullfile(outputDir, sprintf('plumeTrack_output_%s.mat',datestr(now,'YYYY-mm-dd')));
    filename=fullfile(outputDir, 'plumeTrack_output.mat');
    fprintf('Writing parameter file:\n\t%s\n',filename);
    if(exist(filename, 'file')==2)
        delete(filename);
    end
    if(~isempty(Tout))
        T = Tout;
        save(filename,'T','Ref','output_time','-v7.3')
    end
end
fprintf('Done\n');

end

function [listImg,nbFrame,ref,deb,fin]=parseInputDir(inputDir,heads,deb,ref,fin,dN)
%   inputDir = full path to data directory
%   deb = 3 (or 4?) digit integer value of first file to use. Will look for this in the
%        file name, and will use as the reference image. For now lets
%        assume an underscore separates the file index at the end of the
%        filename
%   ref = index of reference image
%   fin = equivalent to deb, OPTIONAL last file to use. Due to shitty index
%        handling, will currently actually grab on more AFTER this
%   dN  = step dN images each iteration
%
% OUT:  listImg = cell containing two columns: {file index, file name},
%                 optionally all other parameter columns if IRBIS ascii table is found
%       nbFrame = number of frames
%

    if nargin<4
        fin = [];
    end

    % Get list and number of frames
    if exist(heads,'file')
        disp('Parameter file found - loading...')
        load(heads); 
        listImg = T;
        
        idxName = cellfun(@str2num,T.Properties.RowNames);
    else
        disp('No parameter file found, assembling file list...')
        listImg = dir(fullfile(inputDir,'*.mat'));
        nbFrame = size(listImg,1);
        
        % >> idxName will correspond to file names, need an equivalent mapping to
        % array indices
        % >> idxNum will correspond to the matlab array indices
        idxName   = zeros(nbFrame,1);
        for ii=1:nbFrame
            idxName(ii) = name2index(listImg(ii).name);
        end
        
        listImg = sortrows([num2cell(idxName) {listImg.name}'],1);
    %     File = sortrows(idxName_cut listImg(inrange)]);
        File = listImg(:,2);
        listImg = table(File,'RowNames',cellstr(string(listImg(:,1))));
        
    end
    idxName = sort(idxName);
    
    % Check for VECTOR deb input
    if ~isscalar(deb)
        disp('Vector frame index input, using list...')
        [idxName_cut,inrange,~]=intersect(idxName,deb);
        fin = max(idxName_cut);
        deb = min(idxName_cut);
%         idxvec = deb;
%         deb = idxvec(1);
%         fin = idxvec(end);
    else
        disp('Scalar frame index input, building list...')
        % Check for empty start/stop indices
        if isempty(deb); deb=min(idxName); end
        if isempty(fin); fin=max(idxName); end

        % Get new index range
        inrange = find(and(idxName>=deb,idxName<=fin));
        inrange = inrange(1:dN:end);
        idxName_cut = idxName(inrange);

        % Check for start/stop not in new range
        if ~ismember(deb,idxName_cut); deb=min(idxName_cut); end
        if ~ismember(fin,idxName_cut); fin=max(idxName_cut); end
        
        
    end
    % Set ref - check empty, check in original range
    if isempty(ref) || ~ismember(ref,idxName); ref=deb; end
    refi    = idxName==ref;
    ref = struct('file',listImg.File{refi},'par',listImg(refi,:));
%     ref.file = listImg.File{refi};
%     ref.par  = listImg(refi,:);

    % Update image list and numframes
    listImg = listImg(inrange,:);    
    nbFrame = size(listImg,1);
end

function idx = name2index(ifile)
%  ifile = full path to file name, including extention
    [~,name,~] = fileparts(ifile);
    istr = split(name,'_');
    idx = double(istr(end));
end

function plotHW(hwind)
% Rows: pt1, pt2.  Columns: WidthX, HeightX, WidthZ, HeightZ


hold on
plot(hwind(5),hwind(6),'wo')
% plot(hwind(:,1:2),hwind(:,3:4),'ro');
% errorbar(x,y,ey,ey,ex,ex,'w.')
end

function plotContourImage(im,outline,cax)
    if nargin<3
        cax=[];
    end
    imagesc(im+max(im(:)).*outline)
    if ~isempty(cax)
        caxis(cax)
    end
    colormap(plasmagrey(200)) % Revert to a classic colormap if you don't have the thermal one
    
end

% function plotContourImage(im,outline,cax)
%     if nargin<3
%         cax=[];
%     end
%     imagesc(im)
%     if ~isempty(cax)
%         caxis(cax)
%     end
%     colormap(thermal(200)) % Revert to a classic colormap if you don't have the thermal one
%     hold on
%     plot(outline(:,2),outline(:,1),'r','LineWidth',1.4)
% end
