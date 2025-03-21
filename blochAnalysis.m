function [k,w,v] = blochAnalysis(M,K,L,pbc,nk,nm)

% define wavenumbers
k = linspace(0,2*pi/L,nk);
nk = length(k);

% initialization
nM = size(K,1);
w = zeros(nm,nk);
v = zeros(nM,nm,nk);

% loop on wave numbers to compute Bloch modes
for i1 = 1:nk
    Q = periodicBoundaryConditions(nM,pbc,k(i1)*L);
    Kq = Q'*K*Q;
    Mq = Q'*M*Q;
    [vv,d] = eigs(Kq,Mq,nm,1e-6);
    v(:,:,i1) = Q*vv;
    w(:,i1) = sqrt(abs(diag(d)));
end
