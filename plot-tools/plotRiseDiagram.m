function varargout = plotRiseDiagram(z,t,D,varargin)
%  R = plotRiseDiagram(z,t,D,varargin)
% IN:
%  z = z vector
%  t = time vector
%  D = thermal data cube struct (assumes plume masks are present)
%  
%
% N/V pairs:
%   mask = data cube mask
%   idx  = subscript selection in 3rd dimension of D
%   tracks = track output from getThermStats
%   detections = output struct from getThermSource to show detected pulses,
%               OR [nx2] vector - first column is trigger time corresponding
%                   to t vector, second column is window height in pixels.
%   atmo = atmospheric profile struct from get/fitAtmoProfile (removed from
%           input data)
%
% R  = output rise diagram image (optional output)


p = inputParser;
addParameter(p,'mask',[])
addParameter(p,'idx',1:size(D,3))
addParameter(p,'cax',[])
addParameter(p,'tracks',[])
addParameter(p,'detections',[])
addParameter(p,'atmo',[])
addParameter(p,'colormap',thermgray(250))
parse(p,varargin{:})

S = p.Results;

assert(length(z)==size(D,1),'z length must match first dimension of D)')
assert(length(t)==size(D,3),'t length must match third dimension of D)')


if isempty(S.idx)
    S.idx = 1:size(D,3);
end
D = D(:,:,S.idx);

if and(~isempty(S.atmo),isstruct(S.atmo))
    D = D - S.atmo.Tinterp - S.atmo.T_halfMax;
end

if ~isempty(S.mask)
    D = D.*S.mask(:,:,S.idx);
end

% if and(~isempty(atmo),isstruct(atmo))
R = squeeze(max(D,[],2));
% else
%     R = squeeze(max(D.T(:,:,S.idx).*D.mask(:,:,S.idx),[],2));
% end
pt = t(S.idx);

if isempty(S.cax)
    S.cax = [min(R(R~=0)) max(R(R~=0))];
end



figure
pcolor(pt,z,R)
shading flat
% view([0 0 1])
colormap(S.colormap)
cb = colorbar;
cb.Label.String = 'Max Temperature [K]';
xlabel('time [s]')
zlabel('z [m]')
caxis(S.cax)
set(gca,'FontSize',12)
axis tight

if ~isempty(S.tracks)
    hold on
    co = get(gca,'ColorOrder');
    co = get(gca,'ColorOrder');
    co = co([1:2 4:7],:); % Because yellow sucks
    for kk=1:length(S.tracks)


        ci = kk-floor(kk/size(co,1))*size(co,1);
        if rem(kk,size(co,1))==0; ci=size(co,1); end
        plotLineError(S.tracks(kk).t,mean(diff(z))*(S.tracks(kk).clustI(:,[1 3])-1)+z(1),co(ci,:)*0.7,0.3,true)
        plot(S.tracks(kk).t, S.tracks(kk).clustZ,'Color',co(ci,:),'LineWidth',1.6)

    end
end
if ~isempty(S.detections)
    hold on
    tS = (t(S.detections.iTrig(:,1))*[1 1])';
    zS = repmat(S.detections.T0.zI([1 end],1),[1 size(S.detections.iTrig,1)]);
    plot(tS,z(zS),'Color',[0   0.447   0.741],'LineWidth',2.5)
end

if nargout==1
    varargout{1} = R;
end
end