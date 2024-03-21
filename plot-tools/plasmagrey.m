function cols = plasmagrey(ncol,greycut)
% Plasma colormap, with a fade to grey at the low end of the colorscale.
% ncol = number of colors
% greycut = fraction [0 - 1], cutoff between grey and colored scales

if nargin<1
    ncol = [];
end
if nargin<2
    greycut = [];
end
if isempty(ncol)
    ncol = 100;
end


    col0 = [0 0 0
            125 125 125
            13 8 135
            84 2 163
            139 10 165
            185 50 137
            219 92 104
            244 136 73
            254 188 43
            240 249 33];

if isempty(greycut)
    cpos = linspace(0,1,length(col0));
else
    assert(and(0<=greycut,greycut<=1),'greycut must be in the range [0,1]')
    cpos = [0 linspace(greycut,1,length(col0)-1)];
end
%     cpos = [0 0.11 0.22 0.33 .44 .55 0.66 .77 0.88 .1];

    cols = plasma(ncol,cpos,col0);

end