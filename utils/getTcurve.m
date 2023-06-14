function [z,T,satfrac,filtfrac] = getTcurve(tk,structName,Tcurve,useZmax) % idx?
% [z,T] = getTcurve(tk,structName,Tcurve)
%                'structName' opts: 'plumeImg'
%                                   'trackImg'
%                                   'trackStats'
%                                   'trackStatsS'
%                                   'trackStatsInt'
%                                   'trackStatsIntS'
%  
%                 'curveName' opts: 'prctile' - 95th percentile
%                                   'mean' 
%                                   'max'

    if nargin<4
        useZmax = false;
    end

    stats = tk.dat.(structName);
%     if isfield(stats,'zmax')
%         z = stats.zmax;
%     else
%         z = stats.z;
%     end

    switch Tcurve
        case 'prctile'
            T = stats.prctile(:,5);
        case 'mean'
            T = stats.mean;
        case 'max'
            T = stats.max;
    end
    
    % Hmmm what if...
%     z = stats.z;
    if isfield(stats,'zmax') && strcmp(Tcurve,'prctile') && useZmax
        z = stats.zmax;
    else
        z = stats.z;
    end
    
    
    satfrac = stats.saturation./stats.N;
    filtfrac = stats.filtFraction;

end