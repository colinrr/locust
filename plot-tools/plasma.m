function cols = plasma(ncol,cpos,col0)

if nargin<1
    ncol = [];
end
if nargin<2
    cpos = [];
end
if nargin<3
    col0 = []; %{'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f','#000000'};
end

if isempty(ncol)
    ncol = 100;
end
if isempty(cpos)
    cpos = [0 0.14 0.29 0.43 0.57 0.71 0.86 1];
end
if isempty(col0)
    col0 = [13 8 135
            84 2 163
            139 10 165
            185 50 137
            219 92 104
            244 136 73
            254 188 43
            240 249 33];
end

    crgb = col0./255; %zeros(size(col0,1),3);
%     for ii=1:length(col0)
%         crgb(ii,:) = hex2rgb(col0{ii});
%     end

    cols = zeros(ncol,3);
    ci = round(cpos.*ncol);
    for ii=1:(length(ci)-1)
        ns = diff(ci(ii:ii+1));
        cols(ci(ii)+1:ci(ii+1),1) = linspace(crgb(ii,1),crgb(ii+1,1),ns)';
        cols(ci(ii)+1:ci(ii+1),2) = linspace(crgb(ii,2),crgb(ii+1,2),ns)';
        cols(ci(ii)+1:ci(ii+1),3) = linspace(crgb(ii,3),crgb(ii+1,3),ns)';
    end
%     cols = flipud(cols);
end