function plotNetwork(n,X,T,def,k,coef)

% constants
if nargin<=4 || isempty(def) 
    def = zeros(size(X,1)*3,1);
    k = 0;
    backstyle = '-k';
else
    backstyle = ':k';
end
if nargin<7 || isempty(coef)
    coef = 1;
end
Nm = size(def,2);
xx = X(unique(T(:)),:);
normref = mean(max(xx,[],1)-min(xx,[],1));

% expand mesh
[X,T,def] = doubleNetwork(n,X,T,[max(X(:,1))-min(X(:,1)) 0],def,k);
Ne = size(T,1);

% loop on modes
for i2 = 1:Nm

    % rescaling the deformed shapes
    u = def(1:3:end,i2);
    v = def(2:3:end,i2);
    mnorm = (norm(u)+norm(v))/2/normref/coef;
    u = u/mnorm;
    v = v/mnorm;
    % theta = def(3:3:end,i2);
    % theta = theta/norm(theta);

    % plot reference configuration
    figure; hold on;
    for i1 = 1:Ne
        Te = T(i1,:);
        plot(X(Te,1),X(Te,2)',backstyle,'linewidth',1)
    end

    % plot deformed configuration
    for i1 = 1:Ne
        Te = T(i1,:);
        plot(X(Te,1)+u(Te),X(Te,2)+v(Te),'-kx','linewidth',2);
    end

end

Lx = max(X(:,1))-min(X(:,1));
Ly = max(X(:,2))-min(X(:,2));
set(gca,'visible','off','PlotBoxAspectRatio',[Lx Ly 1]);