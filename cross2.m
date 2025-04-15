function v = cross2(a,b,dim)
if nargin<2
    dim = 1;
end
if size(a,dim)~=2 || size(b,dim)~=2
    error('input vectors should be 2D')
end
if length(size(a))>2 || length(size(b))>2
    error('only accepts input vectors or matrices')
end
if dim==1
    a = a';
    b = b';
end
v = a(:,1).*b(:,2)-a(:,2).*b(:,1);
if dim==1
    v = v';
end