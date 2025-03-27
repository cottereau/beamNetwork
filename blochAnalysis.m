function [k,w,v] = blochAnalysis(M,K,X,L,pbc,nk,nm)

% define wavenumbers
k = linspace(0,pi/L,nk);
nk = length(k);

% initialization
nM = size(K,1);
w = zeros(nm,nk);
v = zeros(nM,nm,nk);

% loop on wave numbers to compute Bloch modes
tol = 1e-6;
for i1 = 1:nk
    Q = periodicBoundaryConditions(nM,pbc,k(i1)*L);
    Kq = Q'*K*Q;
    Mq = Q'*M*Q;
    [vv,d] = eigs(Kq,Mq,nm,tol);
    d = diag(d);
    v(:,:,i1) = Q*vv;
    w(:,i1) = sqrt(abs(d));
end
% orthonormalization for repeated eigenvalues
Xt = repmat(X(:,1),[1 3])';
v = v./exp(k(i1)*Xt(:));
for i1 = 1:nk
    [v(:,:,i1),~] = qr(v(:,:,i1),0);
end