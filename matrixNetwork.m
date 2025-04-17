function [K,M,Xg,Tg,indT,dKr,dKl,dMr,dMl,indg] = matrixNetwork(varargin)

v = varargin;

switch v{1}
    case 'beam'
        [K,M,Xg,Tg,indT,dKr,dKl,dMr,dMl,indg] = beamNetwork(v{2},v{3},v{4},v{5},v{6},v{7},v{8});

    otherwise 
        error ('unknown link type')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEAM NETWORK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,M,Xg,Tg,indT,dKr,dKl,dMr,dMl,indglob] = beamNetwork(X,T,N,E,rho,S,I)

% constants
nLink = size(T,1);
[nodes,~,indT] = unique(T(:));
indT = reshape(indT,nLink,2);
nNodes = length(nodes);

% full node matrix and connectivity
Xg = [X(nodes,:);zeros(nLink*(N-1),2)];
nXg = size(Xg,1);
Tg = zeros(N*nLink,2);

% allocate matrices
indglob = zeros(nLink,(N+1)*3);
K = sparse(nXg*3,nXg*3);
M = sparse(nXg*3,nXg*3);
dKr = zeros((N+1)*3,(N+1)*3,nLink);
dKl = zeros((N+1)*3,(N+1)*3,nLink);
dMr = zeros((N+1)*3,(N+1)*3,nLink);
dMl = zeros((N+1)*3,(N+1)*3,nLink);

% prepare prompts
blkprompt = repmat(',R',[N 1])';
blkprompt = ['Rg = blkdiag(R' blkprompt(:)' ');'];
dblkprompt = repmat(',dR',[N 1])';
dblkprompt = ['dRg = blkdiag(dR' dblkprompt(:)' ');'];

% loop on links
for i1 = 1:nLink

    % geometry of the link
    X1 = X(T(i1,1),1:2);
    X2 = X(T(i1,2),1:2);
    [theta,L] = cart2pol(X2(1)-X1(1),X2(2)-X1(2));
    Xg((nNodes+(i1-1)*(N-1))+(1:N-1),:) = X1 + ((1:N-1)/N)'*(X2-X1);
    Tg((i1-1)*N+(1:N),:) = [indT(i1,1) (nNodes+(i1-1)*(N-1))+(1:N-1);
                            (nNodes+(i1-1)*(N-1))+(1:N-1) indT(i1,2)]';

    % stiffness and mass matrix of the link in local frame
    [K0,M0,dK0,dM0] = beamLink(N,E,rho,S,I,L);

    % stiffness and mass matrices of the link in global frame
    R = sparse([cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]);
    eval(blkprompt);
    Ki = Rg'*K0*Rg;
    Mi = Rg'*M0*Rg;

    % stiffness and mass perturbations of the link in global frame
    dR = sparse([-sin(theta) cos(theta) 0;-cos(theta) -sin(theta) 0;0 0 0] );
    eval(dblkprompt);
    dKr(:,:,i1) = (Rg'*K0*dRg)+(dRg'*K0*Rg);
    dKl(:,:,i1) = Rg'*dK0*Rg;
    dMr(:,:,i1) = (Rg'*M0*dRg)+(dRg'*M0*Rg);
    dMl(:,:,i1) = Rg'*dM0*Rg;

    % assemble in global matrices
    indloc = [1:3 3*N+(1:3) 3+(1:3*(N-1))];
    indg = [(indT(i1,1)-1)*3+(1:3) ...
               (indT(i1,2)-1)*3+(1:3) ...
               3*(nNodes+(i1-1)*(N-1))+(1:3*(N-1))];
    K(indg,indg) = K(indg,indg) + Ki(indloc,indloc);
    M(indg,indg) = M(indg,indg) + Mi(indloc,indloc);
    indglob(i1,indloc) = indg;

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,M,dK,dM] = beamLink(n,E,rho,S,I,L)
[K0,M0,dK0,dM0] = beamUnit(E,rho,S,I,L/n);
N = 3*(n+1);
K = sparse(N,N);
M = sparse(N,N);
dK = sparse(N,N);
dM = sparse(N,N);
for i1 = 1:n
    ind = (i1-1)*3+(1:6);
    K(ind,ind) = K(ind,ind)+K0;
    M(ind,ind) = M(ind,ind)+M0;
    dK(ind,ind) = dK(ind,ind)+dK0;
    dM(ind,ind) = dM(ind,ind)+dM0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,M,dK,dM] = beamUnit(E,rho,S,I,L)
EIL3 = E*I/L^3;
SL2I = S*L^2/I;
K = EIL3*[ SL2I    0      0 -SL2I    0      0
              0   12    6*L     0  -12   6*L
              0  6*L  4*L^2     0 -6*L 2*L^2
          -SL2I    0      0  SL2I    0      0
              0  -12   -6*L     0   12    -6*L
              0 6*L 2*L^2     0  -6*L  4*L^2];
K = sparse(K);
rSL = rho*S*L/420;
M = rSL*[140    0     0  70     0     0
           0  156  22*L   0    54  -13*L
           0 22*L 4*L^2   0  13*L -3*L^2
          70    0     0 140     0     0
           0   54  13*L   0   156  -22*L
           0 -13*L -3*L^2   0  -22*L 4*L^2 ];
M = sparse(M);
dK = -[1 0 0 1 0 0
       0 3 2 0 3 2
       0 2 1 0 2 1
       1 0 0 1 0 0
       0 3 2 0 3 2
       0 2 1 0 2 1];
dK = sparse(dK).*K;
dM =  [1 0 0 1 0 0
       0 1 2 0 1 2
       0 2 3 0 2 3
       1 0 0 1 0 0
       0 1 2 0 1 2
       0 2 3 0 2 3];
dM = sparse(dM).*M;
end
