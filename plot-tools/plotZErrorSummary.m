function plotZErrorSummary(tt,Z,zErrMax,zErrUncertainty)
% Simple summary plot for zErrorEstimation
        % tt = time or index vector for frames
        % Z  = height vector
        % zErrMax = output direct from zErrorEstimation - zError based on
        %           plume radius dimension
        % zErrUncertainty = Possible range of zErrMax based on uncertainty
        %           in plume centerline position (optional)

        if nargin<4
            zErrUncertainty = squeeze(zErrMax(:,3,:)-zErrMax(:,1,:));
        end
        
        figure('name','Height Error estimation','position',[50 50 800 1000])
        subplot(2,1,1)
        pcolor(tt,Z,squeeze(zErrMax(:,2,:)))
        shading flat
        colormap(copper(200))
        hold on
        elvls = 50:50:(round(max(max(squeeze(zErrMax(:,2,:)))))-50);
        [C,h]=contour(tt,Z,squeeze(zErrMax(:,2,:)),elvls,'w');
        clabel(C,h,'Color','w');
        colorbar
        title('Estimated height error at plume centerline due to projection geometry [m]')
        xlabel('time [s]')
        ylabel('height [m above target]')

        subplot(2,1,2)
        pcolor(tt,Z,zErrUncertainty)
        shading flat
        colormap(copper(200))
        hold on
        elvls = 50:50:(round(max(zErrUncertainty(:)))-50);
        [C,h]=contour(tt,Z,zErrUncertainty,elvls,'w');
        clabel(C,h,'Color','w');
        colorbar
        title('Estimated eight error uncertainty due to plume distance [m]')
        xlabel('time [s]')
        ylabel('height [m above target]')
        
end