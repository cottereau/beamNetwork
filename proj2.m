function p = proj2(u,K,v,ind)
Ne = size(ind,1);
p = zeros(size(u,2));
% p = zeros(1,size(u,2));
for i1 = 1:Ne
    p = p + u(ind(i1,:),:)'*K(:,:,i1)*v(ind(i1,:),:);
%    p = p + sum(conj(u(ind(i1,:),:)).*(K(:,:,i1)*v(ind(i1,:),:)),1);
end
p = p';