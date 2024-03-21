function Poly = designPolygon(img,polyBool,cax)
% imgpath = image array, or cell array of image arrays
            % Designed polygons will be plotted on top of each image
% npolys  = % Boolean vector of length(number of polygons). 
%            true = polygon contains true, false = polygon contains false;
%
    % THE CLICKING BIT
    % Commands: 
    %   left click: add point
    %   right click: delete point
    %   z: zoom 2x towards cursor
    %   x: zoom out 2x away from cursor
    %   y: quit and save poly
    %   q: quit and cancel poly
    %
% C Rowell 2019

if nargin<3
    cax = [];
end
        
    if isnumeric(img) % Single image
        img = {img};
    else
        assert(iscell(img),'Input "img" should be a single image array or cell array of images.')
    end
    nimg = numel(img);
%     np = numel(polyBool);  


    pf = figure('position',[50 200 1800 1000]);
    for mm=1:nimg     
        pax(mm) = tightSubplot(2,2,mm,0,0);
%             if ismember(num2str(despol_idx(mm)),Tout.Properties.RowNames)
%                 imagesc(Frame + Frame.*Tout{num2str(despol_idx(mm)),'Outline'}{1});
%             else
        imagesc(img{mm})
%             end
        axis off
        hold on
        if ~isempty(cax);caxis(cax);end
%         axis equal
        daspect([1 1 1])
    end       
    colormap(thermgray(150))

    pf2=figure('position',[50 200 1500 1000]);
%     load(fullfile(matDir,T.File{num2str(despol_idx(1))}))
%         if ismember(num2str(despol_idx(1)),Tout.Properties.RowNames)
%             imagesc(Frame + Frame.*Tout{num2str(despol_idx(1)),'Outline'}{1});
%         else
    imagesc(img{1})
%         end
    axis off
    hold on
    colormap(thermgray(150))
    if ~isempty(cax);caxis(cax);end
%         linkaxes([pax gca],'xy')
    axis equal   


    polygon = getPoly([pax gca]);
%         if ~npolys(nn)
%             polygon = ~polygon;
%         end
    Poly.P          = polygon;
    Poly.Type       = polyBool;
%     Poly.Design_Idx = despol_idx;
%     Poly.Idx        = IdxPoly{nn};
    Poly.Mask       = poly2mask(polygon(:,1),polygon(:,2),size(img{1},1),size(img{1},2));

    if ~Poly.Type
        Poly.Mask   = ~Poly.Mask;
    end
%         Polys(nn).Outline    = sparse(edge(polygon,'sobel')); 


%     clf(pf)
%     clf(pf2)
%         end

    % Convert poly to mask and add to a cell/struct array

%     save(PolyFile,'Polys')
    close(pf)
    close(pf2)

 %%       
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
    

%         for nn=1:np
%             despol_idx = DesPol_Idx{nn};
%             for mm=1:numel(despol_idx)
%                 load(fullfile(matDir,T.File{num2str(despol_idx(mm))}))      
%                 pax(mm) = tightSubplot(2,2,mm,0,0);
%                 if ismember(num2str(despol_idx(mm)),Tout.Properties.RowNames)
%                     imagesc(Frame + Frame.*Tout{num2str(despol_idx(mm)),'Outline'}{1});
%                 else
%                     imagesc(Frame)
%                 end
%                 axis off
%                 hold on
%                 caxis([190 350])
%         %         axis equal
%                 daspect([1 1 1])
% 
%             end
%             colormap(thermgray(150))           

%             pf2=figure('position',[50 200 1500 1000]);
%             load(fullfile(matDir,T.File{num2str(despol_idx(1))}))
%             if ismember(num2str(despol_idx(1)),Tout.Properties.RowNames)
%                 imagesc(Frame + Frame.*Tout{num2str(despol_idx(1)),'Outline'}{1});
%             else
%                 imagesc(Frame)
%             end
%             axis off
%             hold on
%             colormap(thermgray(150))
%             caxis([190 350])
%     %         linkaxes([pax gca],'xy')
%             axis equal        

%             polygon = getPoly([pax gca],size(Frame));
%     %         if ~npolys(nn)
%     %             polygon = ~polygon;
%     %         end
%             Polys(nn).P          = polygon;
%             Polys(nn).Type       = npolys(nn);
%             Polys(nn).Design_Idx = despol_idx;
%             Polys(nn).Idx        = IdxPoly{nn};
%             Polys(nn).Mask       = poly2mask(polygon(:,1),polygon(:,2),size(Frame,1),size(Frame,2));
% 
%             if ~Polys(nn).Type
%                 Polys(nn).Mask   = ~Polys(nn).Mask;
%             end
%     %         Polys(nn).Outline    = sparse(edge(polygon,'sobel')); 
% 
% 
%             clf(pf)
%             clf(pf2)
% %         end
% 
%         % Convert poly to mask and add to a cell/struct array
% 
%         save(PolyFile,'Polys')
%         close(pf)
%         close(pf2)
end

%%
function Pout = getPoly(h)
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
