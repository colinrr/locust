function varargout=plotLineError(x,y,col,alph,alphmode,lw)
% h=plotLineError(x,y,col,alph,alphmode)
% 
% x = x data
% y = y data
% --> Enter y or x (whichever has the error) as an odd-numbered matrix, 1st
% dimension is the number of samples, 2nd dimension gives
%  [upper_bound(s) measured_val lower_bound(s)]
%
% eg1:  Make a horizontal series with errors in the y axis
%       x  =   [x_1; x_2; ... x_N];
%       y  =   [p90_1 p75_1 p50_1 p25_1 p05_1;
%               p90_2 p75_2 p50_2 p25_2 p05_2;
%                ...
%               p90_N p75_N p50_N p25_N p05_N]
%
%       
% eg2:  Make a vertical series with errors in the x axis
%       x  =   [Val_1+1sigma  Val_1 Val_1-1sigma;
%               Val_2+1sigma  Val_2 Val_2-1sigma;
%                ...
%               Val_N+1sigma  Val_N Val_N-1sigma;]
%       y  =   [y_1; y_2; ... y_N];
%
%
% col      = color of central line value
% alph     = alpha value of error regions
% alphmode = turn on true transparency. In false mode, colours are opaque,
%            but adjusted to rgb values based on background color.
%            Default = false.
%
% C Rowell September 2019

if nargin<6
    lw = 1.5;
end
if nargin<5
    alphmode = [];
end
if nargin<4
    alph = [];
end
if nargin<3
    col = [];
end

if isempty(alphmode)
    alphmode = false;
end
if isempty(alph)
    alph = 0.6;
end
if isempty(col)
    col = [0   0.447   0.741];
end

% assert(isvector(x),'x must needs be a vector, Kemosabe...')
% assert(isvector(y),'y, like x, must needs be a vector, Kemosabe...')

if size(y,2)>1
    plot_y_err = true;
    Err = y;
    X = x;
elseif size(x,2)>1
    plot_y_err = false; % plot x error instead
    Err = x;
    X = y;
end

Nchan = size(Err,2);
Nerr = floor(Nchan/2);
Nmod = mod(Nchan,2);

if numel(alph~=Nerr)
    alph = alph(1).^[Nerr:-1:1]; % Exponential decay of alpha value if only 1 is entered
    
end

% Check for NaNs
IN  = any(isnan(Err),2);
Err = Err(~IN,:);
X   = X(~IN);

for kk=1:Nerr
    upper = Err(:,kk);
    lower = flipud(Err(:,Nchan+1-kk));
    
    if alphmode
        ecol = col;
        ealph = alph(kk); % Effective alpha?
    else
        ecol = rgba2rgb(col,alph(kk),get(gca,'Color'));
        ealph = 1;
    end
    
    if plot_y_err
        h(kk) = fill([X;flipud(X)], [upper; lower], ecol,'linestyle','none','facealpha',ealph);
        
    else
        h(kk) = fill([upper; lower], [X; flipud(X)], ecol,'linestyle','none','facealpha',ealph);
    end
    
    hold on
    
end

% Plotting median
% May consider this part optional...
if Nmod==1
    if plot_y_err
        h(end+1) = plot(X,Err(:,Nerr+1),'LineWidth',lw,'Color',col);
    else
        h(end+1) = plot(Err(:,Nerr+1),X,'LineWidth',lw,'Color',col);
    end
end

if nargout==1
    varargout{1} = h;
end

end

function rgb = rgba2rgb(rgb,alpha,rgb_bg)
%
%

if nargin<3
    rgb_bg = [1 1 1];
end

r = rgb(1);
g = rgb(2);
b = rgb(3);

r2 = (1 - alpha) * rgb_bg(1) + alpha * r;
g2 = (1 - alpha) * rgb_bg(2) + alpha * g;
b2 = (1 - alpha) * rgb_bg(3) + alpha * b;

rgb = [r2 g2 b2];



end