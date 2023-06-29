function [Vtrack] = postProcessTracks(Vtrack,trackPar,smoothLength,oFile)
% [Vtrack] = postProcessTracks(Vtrack)
% Run post-processing on a set of tracked features. Performs 2 main steps:
%       1) Smooth track masks. Uses by default a 3x3x3 box smoother on the
%       3 dimensional masks for each track. This ensures smooth statistical
%       tracking and eliminates outliers and sudden target jumps.
%
%       2) Mutually excludes pixels between tracks, such that no two tracks
%       can contain the same pixels in the same frame. Preference is given
%       to features following later in time, since feature fronts are often
%       best resolved.
%
% IN:   Vtrack       = struct output from pulseTrack.
%                       Updates fields:'clustIdx', 'clustI', 'npx'
%       smoothLength = length of cluster mask smoother in the third (time)
%                       dimension
%       oFile        = OPTIONAL path to output file
%
% C Rowell, Sep 2020
%
fprintf('\n========= pulseTrack Post-processing =========\n')

postProcessSmooth = true;
postProcessClip   = true;

if nargin<4
    oFile = [];
end
if ~isempty(oFile) && ischar(oFile)
    writeMode = true;
else
    writeMode = false;
end

if ischar(Vtrack)
    if exist(Vtrack,'file')
        load(Vtrack)
    end
end

    if postProcessSmooth
        disp('Smoothing clusters...')

        for tk=1:length(Vtrack)
            cmask = false([Vtrack(tk).cubeDims([1 2]) Vtrack(tk).N]);
            for ti=1:Vtrack(tk).N
                cm = cmask(:,:,ti);
                cm(Vtrack(tk).clustIdx{ti}) = true;
                cmask(:,:,ti) = cm;
            end
            cmask = smooth3(double(cmask),'box',[3 3 smoothLength])>0.5;
            for ti=1:Vtrack(tk).N
                Vtrack(tk).clustIdx{ti} = find(cmask(:,:,ti));
            end
        end
    end

    if postProcessClip
        disp('Clipping clusters...')
        VtrackPre = Vtrack;
        nT = length(Vtrack);

        for tk=1:nT-1

            for ti=1:Vtrack(tk).N
                for ci = tk+1:nT
                    [tIcheck,tIloc] = ismember(Vtrack(tk).tI(ti),Vtrack(ci).tI);
                    if tIcheck
                        Vtrack(tk).clustIdx{ti} = setdiff(Vtrack(tk).clustIdx{ti},Vtrack(ci).clustIdx{tIloc});
                    end
                end
                [subI,~] = ind2sub(Vtrack(tk).cubeDims(1),Vtrack(tk).clustIdx{ti});
                Vtrack(tk).clustI(ti,:) = [min(subI) mean(subI) max(subI)];
                Vtrack(tk).npx(ti) = numel(Vtrack(tk).clustIdx{ti});
            end
        end

    end

    if writeMode
        fprintf('Writing processed pulseTrack file:\n\t%s\n',oFile)
        save(oFile,'Vtrack','trackPar')
    end
end