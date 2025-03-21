function [K,M,Xg,Tg,indT] = matrixNetwork(varargin)

v = varargin;

switch v{1}
    case 'beam'
        [K,M,Xg,Tg,indT] = beamNetwork(v{2},v{3},v{4},v{5},v{6},v{7},v{8});

    otherwise 
        error ('unknown link type')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEAM NETWORK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,M,Xg,Tg,indT] = beamNetwork(X,T,N,E,rho,S,I)

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
K = sparse(nXg*3,nXg*3);
M = sparse(nXg*3,nXg*3);

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
    [Ki,Mi] = beamLink(N,E,rho,S,I,L);

    % stiffness and mass matrix of the link in global frame
    R = speye(3);
    R(1:2,1:2) = sparse([ cos(theta) sin(theta); -sin(theta) cos(theta)] );
    blkprompt = 'R';
    for i2=1:N; blkprompt = [blkprompt ',R']; end
    blkprompt = ['Rg = blkdiag(' blkprompt ');'];
    eval(blkprompt);
    Ki = Rg'*Ki*Rg;
    Mi = Rg'*Mi*Rg;

    % assemble in global matrices
    indloc = [1:3 3*N+(1:3) 3+(1:3*(N-1))];
    indglob = [(indT(i1,1)-1)*3+(1:3) ...
               (indT(i1,2)-1)*3+(1:3) ...
               3*(nNodes+(i1-1)*(N-1))+(1:3*(N-1))];
    K(indglob,indglob) = K(indglob,indglob) + Ki(indloc,indloc);
    M(indglob,indglob) = M(indglob,indglob) + Mi(indloc,indloc);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,M] = beamLink(n,E,rho,S,I,L)
[K0,M0] = beamUnit(E,rho,S,I,L/n);
N = 3*(n+1);
K = sparse(N,N);
M = sparse(N,N);
for i1 = 1:n
    ind = (i1-1)*3+(1:6);
    K(ind,ind) = K(ind,ind)+K0;
    M(ind,ind) = M(ind,ind)+M0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,M] = beamUnit(E,rho,S,I,L)
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
end
