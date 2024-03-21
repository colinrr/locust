function plotTimeHistogram(H,ITW,t)
% H = output structure from getTimeHistogram
% OPTIONAL:
%    ITW = threshold control points, [3xn]= [Idx; T0; dT]
if nargin<2
    ITW = [];
end
if nargin<3
    t = [];
end

if ~isempty(t)
    assert(all(size(t)==size(H.Idx)),'Time vector must match frame histogram vector.')
    plot_time = true;
    xlab = 'Time (s)';
    x = t;
else
    plot_time = false;
    xlab = 'Frame Index';
    x = H.Idx;
end

fig = gcf;
pcolor(x,H.Tbins,log10(H.Counts));
shading flat
colormap(gray)
xlabel(xlab)
ylabel('Brightness Temperature (K)')

cb = colorbar;
cb.Label.String = 'log_{10}(Counts)';

if ~isempty(ITW)
    if size(ITW,2)>1
        T0v = interp1(ITW(1,:),ITW(2,:),x,'pchip');
        dTv = interp1(ITW(1,:),ITW(3,:),x,'pchip');
    else
        T0v = ITW(2)*ones(size(x));
        dTv = ITW(3)*ones(size(x));
    end
    hold on
    aa=plot(ITW(1,:),ITW(2,:),'o');
    ab=plot(x,T0v);
    ac=plotLineError(x,[T0v-dTv, T0v, T0v+dTv],ab.Color,0.25,true);
    legs = [aa ab ac(2)];
    labs = {'Filter control points', 'Heaviside temp threshold', 'Heaviside step range'};
end
legend(legs,labs)
end