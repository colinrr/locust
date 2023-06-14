function [sh,eh] = plotScaleError(x,y,gof,qc,ax,sym,col,sms,elw)

if nargin<9 || isempty(elw)
    % sms = 70; % Marker Sizes for scatter
    elw = 1.0; % Marker Sizes for scatter
end
if nargin<8 || isempty(sms)
    % sms = 70; % Marker Sizes for scatter
    sms = 50; % Marker Sizes for scatter
end
if nargin<7 || isempty(col)
    col = 'k';
end
if nargin<6 || isempty(sym)
    sym = 'o';
end
if nargin<5 || isempty(ax)
    ax = gca;
end
if nargin<4 || isempty(qc)
    qc = ones(size(x));
end

% Screen plot
% sms = 70; % Marker Sizes for scatter
% elc = 'k'; % errorbar line color
% elw = 1.8;  %errorbar line width

% Print plot
elc = 'k'; % errorbar line color
% elw = 1.0;  %errorbar line width
% dxs = 0.1; % Plot shift for extra runs

if size(x,2)==3
    dx = diff(x,[],2);
    x = x(:,2);
    useDx = true;
else 
    useDx = false;
    dx = [];
end

if size(y,2)==3
    dy = diff(y,[],2);
    y  = y(:,2);
    useDy = true;
else
    useDy = false;
end

axes(ax)
% eh = errorbar(x,y,dy(:,1),dy(:,2),sym,'LineWidth',elw,'MarkerEdgeColor',elc,'Color',elc);
if useDx && useDy
    eh = errorbar(x,y,dy(:,1),dy(:,2),dx(:,1),dx(:,2),sym,'LineWidth',elw,'Color',col,'MarkerEdgeColor','none');
elseif useDy
    eh = errorbar(x,y,dy(:,1),dy(:,2),sym,'LineWidth',elw,'Color',col,'MarkerEdgeColor','none');
elseif useDx
    eh = errorbar(x,y,[],[],dx(:,1),dx(:,2),sym,'LineWidth',elw,'Color',col,'MarkerEdgeColor','none');
else eh = [];
end
sh = scatter(x,y,sms,gof,sym,'filled','MarkerEdgeColor',col);

scatter(x(qc==0),y(qc==0),sms*3,[0.5 0.4 0.4],'o','LineWidth',elw);  % Highlight bad qc picks
% scatter(x(qc==2),y(qc==2),sms*3,[0 0 0],'o','LineWidth',elw);    % Highlight good qc picks

end