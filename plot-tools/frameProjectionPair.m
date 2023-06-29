function varargout = frameProjectionPair(f1paths,geom1,f2paths,geom2,idx,plotmode)
% ax = frameProjectionPair(f1paths,geom1,f2paths,geom2,idx,plotmode)
% Plot two frames with spatial coordinates either side by side or on top of
% each other to compare projections or different images with same projection.
% IN:   f1paths  = cell of length 2 where:
%           f1paths{1} = directory of frames for first image
%           f1paths{2} = path to header table for first image. 
%       geom1 = geometry file for projection, if applicable. Leave empty if
%           spatial vectors gx,gz are already in the frame file (eg. after
%           running 'interpThermal'
%       geom2 = same as geom1 for second image
%       idx     = frame index (table row name, integer)
%       plotmode = 'side': side by side comparison OR
%                  'stack': plot contours and frame outline of second on top
%                           of first
%                   'mesh': same as 'stack', but also overlays second image
%                           pixel mesh
%
% OUT: ax = plot axis/axes
%
% Bear in mind, this currenlty just gives the broad visual comparison. 
% It does not account for the fact that 'pcolor' plots data points at mesh 
% vertices, rather than cell centers. So plotted pixels/boundaries don't 
% line up perfectly with the actual image pixel boundaries.
%
% C Rowell, March 2018
%

narginchk(4,6)
nargoutchk(0,1)

if nargin<6
    plotmode = 'stack';
end
if nargin<5
    idx = [];
end

% Get first image
load(f1paths{2})
T1 = T;
if isempty(idx)
    idx = str2double(T1.Properties.RowNames{1});
end

load(fullfile(f1paths{1},T1.File{num2str(idx)}))
F1 = Frame;
if all([exist('gx','var'),exist('gz','var')])
    [x1,z1] = meshgrid(gx,gz);
    clear gx gz
elseif ~isempty(geom1)
    if ischar(geom1)
        load(geom1);
    elseif isstruct(geom1)
        geom = geom1;
    end
    [px,pz] = meshgrid(1:size(F1,2),1:size(F1,1));
    [x1,z1] = px2m(px,pz,geom);
else
    error('Cannot create mesh coordinates: either missing spatial vector or geometry file.')
end

% Get second image
load(f2paths{2})
T2 = T;

load(fullfile(f2paths{1},T2.File{num2str(idx)}))
F2 = Frame;
if all([exist('gx','var'),exist('gz','var')])
    [x2,z2] = meshgrid(gx,gz);
    clear gx gz
elseif ~isempty(geom2)
    if ischar(geom2)
        load(geom2);
    elseif isstruct(geom2)
        geom = geom2;
    end
    [px,pz] = meshgrid(1:size(F1,2),1:size(F1,1));
    [x2,z2] = px2m(px,pz,geom);
else
    error('Cannot create mesh coordinates: either missing spatial vector or geometry file.')
end

%% Plot image

switch plotmode
    case {'stack','mesh'}
        figure('position',[50 50 50+size(F1,2) 50+size(F1,1)])
        pcolor(x1,z1,F1)
        shading flat
        colormap(thermal(200))
        
        hold on
        
        contour(x2,z2,F2,5,'Color',[0 0.4 0.7])
        p2(1)=plot([x2(1,1) x2(1,end) x2(end,end) x2(end,1) x2(1,1)],...
            [z2(1,1) z2(1,end) z2(end,end) z2(end,1) z2(1,1)],'b','LineWidth',2);
        pl2{1} = 'Frame 2 contour/border';
        
        if strcmp(plotmode,'mesh')
            p2(2) = mesh(x2,z2,x2*0+mean(F1(:)),'EdgeColor',[0.7 0.7 0.7],'FaceColor','None');
            pl2{2} = 'Frame 2 pixel mesh'; 
        end
        
        
        axis equal tight
        xlabel('X [m from center]')
        ylabel('Z [m above camera]')
        legend(p2,pl2,'location','northwest')
        ax = gca;
    case 'side'
        figure('position',[50 50 50+size(F1,2)*2 50+size(F1,1)])
        tightSubplot(1,2,1)
        pcolor(x1,z1,F1)
        shading flat
        colormap(thermal(200))
        
%          if strcmp(plotmode,'mesh')
%             p2(2) = mesh(x2,z2,x2*0+mean(F1(:)),'EdgeColor',[0.7 0.7 0.7],'FaceColor','None');
%             pl2{2} = 'Frame 2 pixel mesh'; 
%         end
%         axis equal tight
        xlabel('X [m from center]')
        ylabel('Z [m above camera]')
        title('Frame 1')
        daspect([1 1 1])
        ax(1) = gca;
        
        tightSubplot(1,2,2)
        pcolor(x2,z2,F2)
        shading flat
        colormap(thermal(200))
        
%          if strcmp(plotmode,'mesh')
%             p2(2) = mesh(x2,z2,x2*0+mean(F1(:)),'EdgeColor',[0.7 0.7 0.7],'FaceColor','None');
%             pl2{2} = 'Frame 2 pixel mesh'; 
%         end
%         axis equal tight
        xlabel('X [m from center]')
        title('Frame 2')
        ax(2) = gca;  
        daspect([1 1 1])
        linkaxes(ax)
    otherwise
        error('Plotmode not recognized.')
end

%% out
if nargout==1
    varargout{1} = ax;
end
end