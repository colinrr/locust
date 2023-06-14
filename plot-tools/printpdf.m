function printpdf(filename,pdir,dimensions,units,dpi)
%       pdir = printpdf(filename,pdir,dimensions,units,dpi)
%   This is a quick funtion to save pdf images for figures
%   INPUT:      filename   = just that, minus extension
%               dimensions = image size, [x y]
%               path       = save directory
%               units      = units for dimensions: 'centimeters' (default)
%                                                  'inches'
%                                                  'points' (1/72")
%                                                  'normalized'
%               dpi        = dots per inch resolution [def 350]
%
% C Rowell, 2012

% To do ---
%  - varargin setup?
%  - 1 single output path is prob best? A bit funky with extension, so w/e
%  - allow [] inputs for defaulting - DONE
%

%% Parsing
narginchk(3,5)

if nargin<3
    dimensions = [];
end
if nargin<4
    units = [];
end
if nargin<5
    dpi = 350;
end
if isempty(units)
    units = 'centimeters';
end   

if isempty(dimensions)
    dimensions = get(gcf,'PaperPosition');
    dimensions = dimensions(3:4);
end
if isempty(units)
    units = 'centimeters';
end

%%

dpi = ['-r' num2str(dpi)];

set(gcf,'paperunits',units,'paperposition',[0 0 dimensions],'paperpositionmode','manual')
set(gcf,'papersize',dimensions)
fprintf('Writing figure: %s\n',fullfile(pdir, [filename '.pdf']))
print('-dpdf',dpi,fullfile(pdir, [filename '.pdf']));

% export_fig(fullfile(pdir, [filename '.pdf']), '-pdf', dpi,'-transparent')