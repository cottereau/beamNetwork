% description of the unit cell
alpha = pi/4;
Xref = [0 0; 1 sin(alpha);2 0]*15e-3;
Tref = [1 2; 2 3]; 

% physical and geometrical parameters of the links
E = 15.3e6;
rho = 1135;
r = 0.75e-3;
S = pi*r^2;
Im = pi*r^4/4;
cp = sqrt(E/rho);
cb = sqrt(E*Im/rho/S);
Fmax = 5000; % max frequency for transmission analysis

% discretization parameters
nRep = 4;  % number of doubling of the unit cell
nk = 30;   % number of wave numbers
n = 50;   % number of elements for each link
nm = 30;   % number of modes to be computed at each wavenumber
nwf = 500; % number of frequencies for the transmission analysis

% repeating the unit cell
Lx = max(Xref(Tref(:),1))-min(Xref(Tref(:),1));
[X,T] = doubleNetwork(nRep,Xref,Tref,[max(Xref(:,1))-min(Xref(:,1)) 0]);

% boundary conditions
indLeft = Xref(:,1)==min(Xref(:,1));
indLeftRight = (Xref(:,1)==min(Xref(:,1))) | (Xref(:,1)==max(Xref(:,1)));
[~,~,ic] = unique(Tref(:));
Tcounts = accumarray(ic,1);
indCross = Tcounts>2 | (Tcounts>1&indLeftRight);
indMove = indCross&~indLeftRight;
pbc = leftRightPairs(Xref);

% % perturbing the network
% dX = 0.2*2*(rand(size(Xref))-1/2)*15e-3;
% X = Xref;
% X(indMove,:) = X(indMove,:)+dX(indMove,:);
% plotNetwork(0,X,T)
% L0 = vecnorm(Xref(T(:,1),:)-Xref(T(:,2),:),2,2);
% L = vecnorm(X(T(:,1),:)-X(T(:,2),:),2,2);
% dPsi = (L-L0)./L0;
% 
% Bloch analysis of the unperturbed cell
[Kref,Mref,Xgref] = matrixNetwork('beam',Xref,Tref,n,E,rho,S,Im);
[k,wref,vref] = blochAnalysis(Mref,Kref,Xgref,Lx,pbc,nk,nm);
bgref = plotDispersionCurveNetwork(wref,vref,Mref,[],k,cp,cb);
%plotNetwork(0,Xref,T,real(vref(:,2,4)),k(10))

% transmission analysis in frequency for the full length model
[K,M,Xg] = matrixNetwork('beam',X,T,n,E,rho,S,Im);
% [Ut,wt] = transmissionAnalysis(K,M,Xg,Lx,Fmax,nwf);
% % mark in grey shade the potential band gaps
% figure; semilogy(wt/2/pi,abs(Ut))
% addBandGaps(bgref,Ut)

% eigenvalue analysis of the full length model
[v,d] = eigs(K,M,1000,1);
w = sqrt(abs(diag(d)));
hold on; scatter(pi/Lx,w/2/pi,'k')
set(gca,'ylim',[0 5000],'xlim',[0 pi/Lx])

