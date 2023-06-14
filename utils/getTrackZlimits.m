function [RzLims,TzLims,Rqc,Tqc] = getTrackZlimits(trackSet,Rcurve,Tstruct,useZmax)
% Load track z limits for scaling
%   trackSet = list of n tracks for which to obtain z limits
%   
%   RzLims   = [n x 2] array of upper and lower bounds to use for radius data
%   TzLims   = [n x 2] array of upper and lower bounds to use for Temperature data

    xlfile  = '~/Kahuna/data/sabancaya_5_2018/pulseTrack_analysis/pulseTrack_processingV3.xlsx';
    xT = readtable(xlfile,'Sheet','Z_limits_current','ReadVariableNames',true,'ReadRowNames',true);
    
    % Parse sI fields
    if isempty(trackSet)
        trackSet = 1:size(xT,1);
    end
    
    % R Z limits
    if contains(Rcurve,'npx')
        if ~strcmp(Rcurve,'npxInt')
            warning('R Z limits were chosen for interpolated Z values. There may be some issues if using for other npx curves.')
        end
        RzLims = [xT.Rz_npxInt_1(trackSet) xT.Rz_npxInt_2(trackSet)];
        
    elseif contains(Rcurve,'combined') || contains(Rcurve,'Combined')
        RzLims = [xT.Rz_combined_1(trackSet) xT.Rz_combined_2(trackSet)];
        
    else
        error('R z limits are implemented only for npx and combined curves')
    end
    
    % T Z limits
    switch Tstruct
        case 'plumeImg'
            TzLims = [xT.Tz_plumeImg_1(trackSet) xT.Tz_plumeImg_2(trackSet)];
            
        case 'trackImg'
            TzLims = [xT.Tz_trackImg_1(trackSet) xT.Tz_trackImg_2(trackSet)];
            
        case {'trackStats','trackStatsS','trackStatsInt','trackStatsIntS'}
            if ~strcmp(Tstruct,'trackStatsInt')
                warning('T Z limits were chosen for interpolated track values (trackStats). There may be some issues if using for other track curves.')
            end
            if useZmax
                TzLims = [xT.Tz_trackInt_zmax_1(trackSet) xT.Tz_trackInt_zmax_2(trackSet)];
            else
                TzLims = [xT.Tz_trackInt_1(trackSet) xT.Tz_trackInt_2(trackSet)];
            end
            
        otherwise
            error('Tstruct not recognized')
    end
    
    Rqc = xT.QC_Rz(trackSet);
    Tqc = xT.QC_Tz(trackSet);
end