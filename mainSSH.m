% description of the unit cell
alpha = 0;
Xref = [0 0; 1 sin(alpha);2 0]*15e-3;
T = [1 2; 2 3]; 

% physical and geometrical parameters of the links
E = 15.3e6;
rho = 1135;
r = 0.75e-3;
S = pi*r^2;
Im = pi*r^4/4;
cp = sqrt(E/rho);
cb = sqrt(E*Im/rho/S);
Fmax = 1500; % max frequency for transmission analysis

% discretization parameters
nRep = 0;  % number of doubling of the unit cell
nk = 30;  % number of wave numbers
n = 5;    % number of elements for each link
nm = 18;  % number of modes to be computed at each wavenumber
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
% I need to construct element-by-element scalar product in order to use
% simply the perturbation formula
[Kref,Mref,Xgref] = matrixNetwork('beam',Xref,T,n,E,rho,S,Im);
[k,wref,vref] = blochAnalysis(Mref,Kref,Xgref,Lx,pbc,nk,nm);
bgref = plotDispersionCurveNetwork(wref,vref,[],k,cp,cb);
%plotNetwork(0,Xref,T,real(vref(:,2,4)),k(10))

% transmission analysis in frequency for the unperturbed cell
[Utref,wref] = transmissionAnalysis(Kref,Mref,Xgref,Lx,Fmax,nwf);
% mark in grey shade the potential band gaps
figure; semilogy(wref/2/pi,abs(Utref))
addBandGaps(bgref,Utref)

% % Bloch analysis of the perturbed cell
% [K,M,Xg,Tg,indT] = matrixNetwork('beam',X,T,n,E,rho,S,Im);
% [k,w0,v0] = blochAnalysis(M,K,Lx,pbc,nk,nm);
% bg = plotDispersionCurveNetwork(w0,v0,[],k,cp,cb);
% 
% % transmission analysis in frequency for the perturbed cell
% [Ut,wt] = transmissionAnalysis(K,M,Xgref,Lx,Fmax,nwf);
% figure;semilogy(wt/2/pi,abs(Ut))
% addBandGaps(bg,Ut)
