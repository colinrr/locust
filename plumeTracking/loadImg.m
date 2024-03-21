function [im,ts] = loadImg( iname, rotateFlag, angle )
%LOADIMG Summary of this function goes here
%   iname       = input file path/name.ext
%   rotateFlag  = [opt] in case images need rotation?
%   angle       = [opt] angle to rotate
% OUT:
%   im          = image array
%   ts          = timestamp

if nargin<2
    rotateFlag=0;
    angle=-90;
else
    if nargin<3
    angle=-90;
    end
end


[~,fname,ex] = fileparts(iname);
im=load(iname);
ind=fieldnames(im) ;
% Assumes an order to the field names, but w/e
ts=getfield(im,'File_DateTime'); % Timestamp
im=getfield(im,'Frame'); % Frame

% Edit timestamp to get s + ms
% ts = [ts(1:5) ts(6)+ts(7)/1000];

% if nargout>1
%     nbFrame=size(dir([name '*.mat']),1);
% end

if rotateFlag
    im=imrotate(im,angle);
end

end