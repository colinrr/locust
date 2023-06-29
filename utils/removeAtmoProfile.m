function [T,mask,threshT] = removeAtmoProfile(mask,T,atmo,nHalfMax,removeBackground)
%   [T,mask] = removeAtmoProfile(mask,T,atmo,nHalfMax,removeBackground)
%   mask cube to trim
%   T        = thermal data set (raw)?
%   atmo     = atmo profile struct
%   nHalfMax = 2, typically. After removing atmospheric profile, a data
%                threshold will be set,
%                T_THRESH = -nHalfMax*(Tmax-T_halfMax), where (Tmax-T_halfMax) is
%                the half-width of the temperature distribution for
%                temperature values that fit the atmospheric profile.
%                (see fitAtmoProfile to make sense of this).
%                Pixels below this value will be discarded from the mask as
%                background.
%   removeBackground = OPTIONAL FLAG (default = TRUE). When true, will set
%                 all T values below T_THRESH EQUAL TO T_THRESH to provide
%                 a uniform temparature background.
%   

    if nargin<4
        nHalfMax = [];
    end
    if nargin<5
        removeBackground = true;
    end
    if isempty(nHalfMax)
        nHalfMax =2;
    end

    threshT = -nHalfMax*(atmo.Tmode-atmo.T_halfMax);

    T = T-atmo.Tinterp-atmo.T_halfMax;
    threshI = T<threshT;
    if removeBackground
        T(threshI) = threshT;
        mask(threshI) = 0;

    end
%     mask(threshI) = 0;


end