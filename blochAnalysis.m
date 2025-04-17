function [k,w,v] = blochAnalysis(M,K,L,pbc,nk,nm)

% define wavenumbers
k = linspace(0,pi/L,nk);
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
    % if i1==1
    %     Kq = Kq(3:end,3:end);
    %     Mq = Mq(3:end,3:end);
    %     Nq = size(Kq,1);
    %     [vv,dd] = eigs(Kq,Mq,nm-2,'smallestabs');
    %     d = [0;0;diag(dd)];
    %     v1 = zeros(Nq+2,1); v1(1:3:end,1) = 1;
    %     v2 = zeros(Nq+2,1); v2(2:3:end,1) = 1;
    %     vv = [v1 v2 [zeros(2,nm-2);vv]];
    % else
        [vv,d] = eigs(Kq,Mq,nm,'smallestabs',tolerance=1e-6);
        d = diag(d);
%    end
    vv = Q*vv;
    v(:,:,i1) = vv*diag(1./sqrt(abs(diag(vv'*M*vv))));
    w(:,i1) = sqrt(abs(d));
end

% % orthonormalization for repeated eigenvalues
% Xt = repmat(X(:,1),[1 3])';
% v = v./exp(k(i1)*Xt(:));
% for i1 = 1:nk
%     [v(:,:,i1),~] = qr(v(:,:,i1),0);
% end