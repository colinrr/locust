function varargout = plotImageMap(gm,img,imgMode)
% ax = plotImageMap(geom,img,mode)
% Show distortion of images projected with 'mapPixels' and 'px2m'.
%   IN:  gm   = geometry structure from mapPixels.m
%   OPTIONAL:
%        img  = path to a single thermal frame image to project
%        mode = 'full' [default] or 'wire'.
%
% OPTIONAL OUTPUT:
%   ax = axes handles for output plots
%
% C Rowell, March 2020

narginchk(1,3)
if nargin<3
    imgMode = 'full';
end
if nargin<2
    img = [];
end

dx      = 0.11;
ppads   = [0.08 0.03 0.09 0.03];
lw      = 1.5;

% Get coordinate vectors and meshes.
[xpg,ypg] = meshgrid(1:gm.im_size(2),1:gm.im_size(1)); % Pixel grid
[xng,yng] = meshgrid(0.5:1:(gm.im_size(2)+0.5),0.5:1:(gm.im_size(1)+0.5));  % Node grid

px = [xpg(:,1); xpg(end,2:end-1)'; flipud(xpg(:,end)); fliplr(xpg(1,2:end-1))'  ];
py = [ypg(:,1); ypg(end,2:end-1)'; flipud(ypg(:,end)); fliplr(ypg(1,2:end-1))'  ];

% Outer boundary coords (nodes)
pxb = [xng(:,1); xng(end,2:end-1)'; flipud(xng(:,end)); fliplr(xng(1,2:end-1))'  ];
pyb = [yng(:,1); yng(end,2:end-1)'; flipud(yng(:,end)); fliplr(yng(1,2:end-1))'  ];

[xc_b,yc_b] = px2m(px,py,gm);
[xn_b,yn_b] = px2m(pxb,pyb,gm);
[xg,yg] = px2m(xng,yng,gm);
[xp,yp] = px2m(xpg,ypg,gm);

%% No Image: Plot pre- and post-projection pixel meshes and frames
if exist(img,'file')
    load(img)
    imgMode = 'full';
elseif all(size(img)==gm.im_size)
    Frame = img;
    imgMode = 'full';
else
    imgMode = 'wire';
end


    switch imgMode
        case 'full'
            figure('position',[100 100 1000 500])
            ax(1) = tightSubplot(1,2,1,dx,[],ppads);
            pcolor(xpg,ypg,Frame)
            shading flat
            title('Original frame (pixels)','Interpreter','Latex')
            set(gca,'YDir','reverse')
%             axis equal
            grid on
            hold on
            xlabel('$i$ (pixels)','Interpreter','Latex')
            ylabel('$j$ (pixels)','Interpreter','Latex')
            daspect([1 1 1])
            
            ax(2) = tightSubplot(1,2,2,dx,[],ppads);
            pcolor(xp,yp,Frame)
            shading flat
            hold on
            plot(xn_b,yn_b,'r')
            title('Projected frame (m)','Interpreter','Latex')
            xlabel('Meters from center','Interpreter','Latex')
            ylabel('Meters above camera','Interpreter','Latex')
            daspect([1 1 1])
            grid on
        case 'wire'

            figure('position',[100 100 1000 500])
            % subplot(1,2,1)
                mesh(xg,yg,xng*0)
                view([0 0 1])
                hold on
            %     scatter3(xpg(:),ypg(:),xpg(:)*0,20,'g.')
                title('pix')
            %     set(gca,'YDir','reverse')
                axis equal
                grid on
                plot3(xn_b,yn_b,xn_b*0,'r')
            % subplot(1,2,2)
            % view([0 0 1])
            % title('m')
            % axis equal
            % grid on

            figure('position',[100 100 1000 500])
            ax(1) = subplot(1,2,1);
            plot(px,py,'.-')
            title('Original image (pixels)')
            set(gca,'YDir','reverse')
            axis equal
            grid on
            hold on
            ax(2) = subplot(1,2,2);
            plot(xc_b,yc_b,'.-')
            hold on
            plot(xn_b,yn_b,'r')
            title('Projected Image (m)')
            xlabel('Meters from center')
            ylabel('Meters above camera')
            axis equal
            grid on
    end

    lp1 = [];
    lp2 = [];
    ll = {};
    if isfield(gm,'ref_pix_ij')
        lp1 = [lp1 plot(ax(1),gm.ref_pix_ij(2),gm.ref_pix_ij(1),'xc','LineWidth',1.5)];
        [xR,zR] = px2m(gm.ref_pix_ij(2),gm.ref_pix_ij(1),gm);
        lp2 = [lp2 plot(ax(2),xR,zR,'xc','LineWidth',1.5)];
        ll = [ll 'Reference pixel'];
    end
    if isfield(gm,'target_pix_ij')
        lp1 = [lp1 plot(ax(1),gm.target_pix_ij(2),gm.target_pix_ij(1),'xw','LineWidth',1.5)];
        [xT,zT] = px2m(gm.target_pix_ij(2),gm.target_pix_ij(1),gm);
        lp2 = [lp2 plot(ax(2),xT,zT,'xw','LineWidth',1.5)];
        ll = [ll 'Projected vent pixel'];
    end
    if ~isempty(ll)
        legend(lp1,ll,'Interpreter','Latex')
        legend(lp2,ll,'Interpreter','Latex')
    end
    
%     if ~isempty(img)
%         warning(sprintf('Could not find image file:\n\t%s',img))
%     end

    if nargout == 1
        varargout{1} = ax;
    end
end