function [X,T,v] = doubleNetwork(n,X,T,d,v,k)

% constants
if nargin<5
    v = zeros(3*size(X,1),1);
    k = 0;
end

% remove unused nodes
[nodes,~,indT] = unique(T(:));
X = X(nodes,:);
T = reshape(indT,size(T));

% loop on number of doubling demanded
for i1 = 1:n
    [X,T,v] = doubleNetworkOnce(X,T,(2^(i1-1))*d,v,k);
end

end
%============================
function [X,T,v] = doubleNetworkOnce(X,T,d,v,k)
% repeat pattern in the chosen direction
T = [T;T+size(X,1)];
X = [X;X+d];
u = [v(1:3:end,:);v(1:3:end,:)*exp(1i*k*d(1))];
v2 = [v(2:3:end,:);v(2:3:end,:)*exp(1i*k*d(1))];
w = [v(3:3:end,:);v(3:3:end,:)*exp(1i*k*d(1))];

% remove repeated nodes
[X,ia,ic] = unique(X,'stable','rows');
T = ic(T);
v = zeros(3*length(ia),size(v,2));
v(1:3:end,:) = u(ia,:);
v(2:3:end,:) = v2(ia,:);
v(3:3:end,:) = w(ia,:);

end