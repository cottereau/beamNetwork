function Q = periodicBoundaryConditions(N,pbc,kx)

% initialization
ind0 = ((pbc(:,1)-1)*3)+(1:3);
ind1 = ((pbc(:,2)-1)*3)+(1:3);

% loop on boundary condition pair
Q = speye(N);

% adding boundary condition along x direction
Q = Q + sparse(ind1(:),ind0(:),exp(1i*kx),N,N);

% removing unnecessary columns
ind = setdiff(1:N,ind1(:));
Q = Q(:,ind);
