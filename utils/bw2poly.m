function B = bw2poly(bw,conn,opts)
% May need testing for cases with holes?

if nargin<2
    conn=8;
end
if nargin<3
    opts = {'holes'};
end

if isempty(conn)
    conn=8;
end

[Bcell,L,n] = bwboundaries(bw,conn,opts{:});

B(length(Bcell)) = struct('Length',[],'X',[],'Y',[],'IsFilled',[]);
for ii=1:length(Bcell)
    B(ii).Length = size(Bcell{ii},1);
    B(ii).X = Bcell{ii}(:,2);
    B(ii).Y = Bcell{ii}(:,1);
    if ii<=n
        B(ii).IsFilled = 1;
    else
        B(ii).IsFilled = 0;
    end
end

end