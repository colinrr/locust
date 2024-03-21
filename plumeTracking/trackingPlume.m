function [ contourOut, closeImg ] = trackingPlume(i, ref, iDir, listImg, dt)
%TRACKINGPLUME Summary of this function goes here
%   INPUT heavily modified by CR...July 2018
%       i       = matlab array index corresponding to starting image in listImg
%       ref     = filename for reference image
%       iDir    = input directory
%       listImg = table containing column vectors: {image index, img file}
%       dt      = step size for images
%
%   OUTPUT: contourOut = sparse matrix giving plume OUTLINE
%           closeImg   = sparse matrix giving plume MASK

% Modified by Colin Rowell, July-Aug 2018
%% Load Images

% CR FIX: set ref image as input, and actual filenames
% as well, so ref image isn't fixed as '1'
% fprintf('Prev: ')
% fprintf('i = %i, Index(%i): %s, Index(%i+1): %s  (dN=%i)\n',i,i,listImg.Properties.RowNames{i},i,listImg.Properties.RowNames{i+1},dt)
src=loadImg(fullfile(iDir,listImg.File{i}));
% fprintf('Next: ')
src2=loadImg(fullfile(iDir,listImg.File{i+1})); % CR: Removed dN and replaced with +1 (input list accounts for dN, now)
% fprintf('Ref:  ')
origin=loadImg(fullfile(iDir,ref));

% TEMP MESS AROUND
origin = ones(size(origin))*min(origin(:));
%% Get a mask of the moving objects
% Load the consecutive differentiation
maskPrec=dif2mask(src,src2);
maskPrec=openAndClean(maskPrec);

% Load the original differentiation
maskOri=dif2mask(origin,src2);
maskOri(maskOri>1)=1;

% Create a mask from both differentiations (moving object which was not
% there at the beginning
mask=maskOri+maskPrec;
mask=imfill(mask,'holes');

%% Clean Image
cleanImg=openAndClean(biggestConnexComponent(mask));

% We redefine the edges since the wavelet transformation is coarse
rIm=imreconstruct(cleanImg,src2);

% We chose the larger area reconstruct inside the mask as the plume.
val=unique(rIm);
m=0;
ind=0;
for j=1:size(val,1)
    m2=size(find(rIm==val(j)),1);
    if m2>m
        m=m2;
        ind=j;
    end
end
val=val(j);

% We binarize the area
rIm(rIm~=val)=-Inf;
rIm(rIm==val)=1;
rIm(rIm~=1)=0;
cleanImg=rIm.*cleanImg;

% We do a closing transformation to remove the last outliers
se=strel('disk',2);
closeImg=imclose(cleanImg,se);
closeImg=imfill(closeImg,'holes');
closeImg=biggestConnexComponent(closeImg);

% We compute the contour to draw it on the original image
BW=edge(closeImg,'sobel'); % Consider using 'sobel' here instead of 'canny'

%%%%%%% COLIN EDITS %%%%%%
% ----- Colin's bad line edit ----
% Get edge as a line
% idx = find(BW); [z,x] = ind2sub(size(BW),idx);
% zo = zeros(size(z)); xo = zeros(size(x));
% zo(1) = z(1); xo(1) = x(1);
% a = plot(xo(1),zo(1),'r.');
% used_idx = zeros(size(z)); used_idx(1) = 1;
% nidx = 1;
% for kk = 1:length(z)-1
%     zn = z(~used_idx);
%     xn = x(~used_idx);
% %     d  = sqrt((xn-x(kk)).^2 + (zn-z(kk)).^2);
%     [~,~,tidx] = closest2d(x(nidx),z(nidx),xn,zn);
%     xo(kk+1) = xn(tidx);
%     zo(kk+1) = zn(tidx);
%     nidx = find(and(xn(tidx)==x,zn(tidx)==z));
%     used_idx(nidx) = 1;
%     delete(a)
%     a = plot(xo(1:kk+1),zo(1:kk+1),'r.');
% end
% [edgelist, ~, ~] = edgelink(BW);
% 
% % Check number of edge lines
% if numel(edgelist)>1
%     warning('More than one edge detected!')
%     oline = edgelist{1};
%     for ii=2:length(edgelist)
%         endpoints = edgelist{ii}([1 end],:);
%         [~,~,ix] = closest2d(oline(end,2),oline(end,1),endpoints(:,2),endpoints(:,1));
%         if ix==1
%             oline = [oline; edgelist{ii}];
%         elseif ix==2
%             oline = [oline; flipud(edgelist{ii})];
%         else
%             error('What did you fuck up?')
%         end
%     end
%     contourOut=oline;
% else
%     contourOut=edgelist{1};
% end
contourOut = sparse(BW);

% -vvv- Original out -vvv-
% contourOut=src2+max(src2(:))*BW;
%%%%% END COLIN EDITS %%%%%%
end
