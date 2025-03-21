function [Utref,wf] = transmissionAnalysis(K,M,X,L,Fmax,nwf)

% define excitation vector
indexc = (find(X(:,1)<=L/20)-1)*3+(1:3); Nexc = nnz(indexc(:));
f = zeros(size(K,1),1); f(indexc(:),:) = rand(Nexc,1);

% define observation vector
indobs = (find(X(:,1)>=L*(19/20))-1)*3+(1:3); Nobs = nnz(indobs(:));
obs = zeros(1,size(K,1)); obs(:,indobs(:)) = rand(1,Nobs);

% define frequency range
wf = 2*pi*linspace(0,Fmax,nwf);

% initialization
Utref = zeros(nwf,1);

% loop on frequency and solution
for i1 = 1:nwf
    Utref(i1) = obs(1,:)*((K-wf(i1).^2*M)\f(:,1));
end
