% description of the unit cell
% square lattice (1D)
% Xref = [0 0; 0 2;2 0;2 2]*15e-3;
% T = [1 2; 1 3;2 4];
% square lattice with one diagonal (1D)
% Xref = [0 0; 0 2;2 0;2 2]*15e-3;
% T = [1 2; 1 3;2 4;1 4]; % one diagonal
% square lattice with two diagonals (1D)
Xref = [0 0; 0 2;2 0;2 2;1 1]*5e-3;
T = [1 2; 1 3;2 4;1 5;2 5;3 5;4 5]; % two diagonals

% physical and geometrical parameters of the links
E = 15.3e6;
rho = 1135;
b = 10e-3;
r = 1.5e-3;
S = b*r;        % square section
Im = b*r^3/12;    % square inertia 
% S = pi*r^2;     % circular section
% Im = pi*r^4/4;  % circular inertia
Fmax = 3000; % max frequency for transmission analysis

% discretization parameters
nRep = 5;  % number of doubling of the unit cell
nk = 30;  % number of wave numbers
n = 5;    % number of elements for each link
nm = 100;  % number of modes to be computed at each wavenumber
nwf = 500;% number of frequencies for the transmission analysis

% repeating the unit cell
[Xref,T] = doubleNetwork(nRep,Xref,T,[max(Xref(:,1))-min(Xref(:,1)) 0]);
Lx = max(Xref(T(:),1))-min(Xref(T(:),1));
Ly = max(Xref(T(:),2))-min(Xref(T(:),2));

% boundary conditions
indLeft = Xref(:,1)==min(Xref(:,1));
indLeftRight = (Xref(:,1)==min(Xref(:,1))) | (Xref(:,1)==max(Xref(:,1)));
[~,~,ic] = unique(T(:));
Tcounts = accumarray(ic,1);
indCross = Tcounts>2 | (Tcounts>1&indLeftRight);
indMove = indCross&~indLeftRight;
pbc = leftRightPairs(Xref);

% Bloch analysis of the unperturbed cell
[Kref,Mref,Xgref,Tgref,~,dKr,dKl,dMr,dMl,indg] = matrixNetwork('beam',Xref,T,n,E,rho,S,Im);
[k,wref,vref] = blochAnalysis(Mref,Kref,Lx,pbc,nk,nm);
bgref = plotDispersionCurveNetwork(wref,vref,Mref,[],k);

% % transmission analysis in frequency for the unperturbed cell
% [Utref,wtref] = transmissionAnalysis(Kref,Mref,Xgref,Lx,Fmax,nwf);
% figure; semilogy(wtref/2/pi,Utref)
% addBandGaps(bgref,Utref)

% perturbing the network
L0 = vecnorm(Xref(T(:,1),:)-Xref(T(:,2),:),2,2);
dX = 0.1*2*(rand(size(Xref))-1/2)*mean(L0(:));
dX(~indMove,:) = 0;
X = Xref+dX;
%plotNetwork(0,X,T)
L = vecnorm(X(T(:,1),:)-X(T(:,2),:),2,2);

% Bloch analysis of the perturbed cell
[K,M,Xg,Tg,indT] = matrixNetwork('beam',X,T,n,E,rho,S,Im);
[~,w0,v0] = blochAnalysis(M,K,Lx,pbc,nk,nm);
bg = plotDispersionCurveNetwork(w0,v0,M,[],k);

% % transmission analysis in frequency for the perturbed cell
% [Ut,wt] = transmissionAnalysis(K,M,Xgref,Lx,Fmax,nwf);
% figure;semilogy(wt/2/pi,Ut)
% addBandGaps(bg,Ut)

% analytical estimation of perturbed dispersion curve
dL = dot(dX(T(:,1),:)-dX(T(:,2),:),Xref(T(:,1),:)-Xref(T(:,2),:),2)./(L0.^2);
dL = reshape(dL,[1 1 length(dL)]);
dphi = cross2(Xref(T(:,1),:)-Xref(T(:,2),:),dX(T(:,1),:)-dX(T(:,2),:),2)./(L0.^2);
dphi = reshape(dphi,[1 1 length(dphi)]);
dK = dKr.*dphi + dKl.*dL;
dM = dMr.*dphi + dMl.*dL;
dw = zeros(nm,nk);
for i1 = 1:nk
    vdKv = real(proj2(vref(:,:,i1),dK,vref(:,:,i1),indg));
    vdMv = real(proj2(vref(:,:,i1),dM,vref(:,:,i1),indg));
    dw(:,i1) = (vdKv-((wref(:,i1).^2).*vdMv))./(2*wref(:,i1));
end
[w0(:,10)-wref(:,10) dw(:,10)]
% for i1 = [1 nk]
%     i1
%     indrep = find(diff(wref(:,i1))<1e-3*mean(diff(wref(:,i1))));
%     v1dKv1 = real(proj2(vref(:,indrep,i1),dK,vref(:,indrep,i1),indg));
%     v1dKv2 = real(proj2(vref(:,indrep,i1),dK,vref(:,indrep+1,i1),indg));
%     v2dKv2 = real(proj2(vref(:,indrep+1,i1),dK,vref(:,indrep+1,i1),indg));
%     v1dMv1 = real(proj2(vref(:,indrep,i1),dM,vref(:,indrep,i1),indg));
%     v1dMv2 = real(proj2(vref(:,indrep,i1),dM,vref(:,indrep+1,i1),indg));
%     v2dMv2 = real(proj2(vref(:,indrep+1,i1),dM,vref(:,indrep+1,i1),indg));
%     theta = 2*(v1dKv2-wref(indrep,i1).^2.*v1dMv2) ...
%                    ./((v1dKv1-v2dKv2)-wref(indrep,i1).^2.*(v1dMv1-v2dMv2));
%     theta = atan(theta);
%     dw(indrep,i1) = (cos(theta).^2.*(v1dKv1-wref(indrep,i1).^2.*v1dMv1) ...
%                     +sin(theta).^2.*(v2dKv2-wref(indrep,i1).^2.*v2dMv2) ...
%        +2*cos(theta).*sin(theta).*(v1dKv2-wref(indrep,i1).^2.*v1dMv2)) ...
%                                                      ./(2*wref(indrep,i1));
%     dw(indrep+1,i1) = (sin(theta).^2.*(v1dKv1-wref(indrep,i1).^2.*v1dMv1) ...
%                     +cos(theta).^2.*(v2dKv2-wref(indrep,i1).^2.*v2dMv2) ...
%        -2*cos(theta).*sin(theta).*(v1dKv2-wref(indrep,i1).^2.*v1dMv2)) ...
%                                                      ./(2*wref(indrep,i1));
% end
plotDispersionCurveNetwork(wref+dw,vref,M,[],k);

