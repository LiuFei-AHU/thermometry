function sz = size(a,dim)

if min(size(a.kmask)) > 1
    sz = [sum(a.kmask(:)) numel(a.kmask)];
else
    sz = [sum(a.kmask)*length(a.kmask) length(a.kmask)*length(a.kmask)]; 
end

if nargin > 1
    sz = sz(dim);
end