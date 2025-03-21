function pbc = leftRightPairs(X)

% constants
indLeft = find(X(:,1)==min(X(:,1)));
yLeft = X(indLeft,2);
indRight = find(X(:,1)==max(X(:,1)));
yRight = X(indRight,2);
N = length(indLeft);

% initialization
pbc = zeros(N,2);

% loop on pairs
for i1 = 1:N
    pbc(i1,:) = [indLeft(i1) indRight(yRight==yLeft(i1))];
end
