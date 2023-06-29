
function Vtrack  = truncateVtrack(Vtrack,N)
% Quick function to crop out pulseTracker results for a subset of frames.
% Vtrack    = output struct from pulseTracker.
% N         = scalar or [1x2]: crops out frame indices outside this range for all
%             fields.

    ff = fieldnames(Vtrack);
    
    for ti=1:length(Vtrack)
    
        for fi = 1:length(ff)
            fn = ff{fi};
            fSz = size(Vtrack(ti).(fn));
            dimChk = fSz==Vtrack(ti).N;
            if any(dimChk)
                for dd=1:ndims(Vtrack(ti).(fn))
                    if dimChk(dd)
                        idx{dd} = 1:N;
                    else
                        idx{dd} = 1:size(Vtrack(ti).(fn),dd);
                    end
                end
    %             if iscell(Vtrack.(fn))
    %                 Vtrack.(fn) = Vtrack.(fn){idx{:}};
    %             elseif isnumeric(Vtrack.(fn))
                Vtrack.(fn) = Vtrack(ti).(fn)(idx{:});
    %             end
            end
        end
        Vtrack(ti).N = N;
    end
end