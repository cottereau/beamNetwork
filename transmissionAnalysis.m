function [Utref,wf] = transmissionAnalysis(K,M,X,L,Fmax,nwf)

% rescale
X = X(:,1)-min(X(:,1));

% define excitation vector
indexc = (find(X(:,1)<=L/20)-1)*3+(1:3); 
Nexc = nnz(indexc(:));
f = zeros(size(K,1),Nexc); 
%f(indexc(:),:) = rand(Nexc,1);
f(indexc(:),:) = eye(Nexc);

% define observation vector
indobs = (find(X(:,1)>=L*(19/20))-1)*3+(1:3); 
Nobs = nnz(indobs(:));
obs = zeros(Nobs,size(K,1));
obs(:,indobs) = eye(Nobs);

% define frequency range
wf = 2*pi*linspace(0,Fmax,nwf)';

% initialization
Utref = zeros(Nobs,Nexc,nwf);

% loop on frequency and solution
for i1 = 1:nwf
    Utref(:,:,i1) = obs*((K-wf(i1).^2*M)\f);
end

% condensate
Utref = squeeze(mean(mean(abs(Utref),1),2));