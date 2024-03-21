function [Fgmask,Poly] = makeForegroundMask(ref,threshFG,polys)
% Simple function to create a mask for image foreground.
% IN: 
%   ref      = path to reference image
%   FGthresh = threshold value. Values BELOW this will be cut
%
% OPTIONAL:
%   polys = use custom polygon add mask elements. Options:
%            To DESIGN polygons: Boolean vector of length(number of polygons). 
%              '->  true  = polygon contains true (adds to mask), 
%              '->  false = polygon contains false (subtracts from mask);
%           'path/img' = path to existing polygon file

narginchk(2,3)
if nargin<3
    polys = [];
end

if islogical(polys)
    design_polys = true;
    load_polys = false;
elseif ~isempty(polys) && ischar(polys)
    design_polys = false;
    load_polys = true;
    load(polys)
else
    design_polys = false;
    load_polys = false;
end

    load(ref)
    Fgmask         = Frame>threshFG;
    cax = [threshFG max(Frame(:))];
    
    if design_polys
        for nn=1:numel(polys)
            Poly(nn) = designPolygon(Frame.*Fgmask,polys(nn),cax);
            
            % Add polygon(s) to mask
            Fgmask = Fgmask.*Poly(nn).Mask;
        end
        
%         disp('hold it...')
    end
    
end