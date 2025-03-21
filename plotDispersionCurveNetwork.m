function bg = plotDispersionCurveNetwork(w,v,AB,kx,cp,cb)

% remove end points
% w = w(:,2:end-1);
% v = v(:,:,2:end-1);
% kx = kx(2:end-1);

% change order of dispersion curves to track modes
for i1 = 1:length(kx)-1
    v1 = v(:,:,i1);
    v2 = v(:,:,i1+1);
    mac = abs(v1'*v2)./sqrt(diag(v1'*v1).*diag(v2'*v2)');
%    figure; surf(mac),shading flat, view(2)
    [~,ind]=max(mac);
    if i1==length(kx)-1
        disp('ok')
    end
    if length(unique(ind))~=size(v,2)
        [ind2,inda] = unique(ind);
        ind = 1:length(ind);
        ind(inda) = ind2;
    end
    v(:,ind,i1+1) = v(:,:,i1+1);
    w(ind,i1+1) = w(:,i1+1);
end

% mark in grey shade the potential band gaps
f = w/2/pi;
minw = min(f,[],2);
maxw = max(f,[],2);
bg = [maxw(1:end-1) minw(2:end)];
for i1 = 1:size(bg,1);bg(2:end,1)=max(bg(1:end-1,1),bg(2:end,1));end
bg = bg(bg(:,2)>bg(:,1),:);

% constants;
nk = size(w,2);
if isempty(AB)   % 1D
    figure; plot(kx,f,'-x','linewidth',1.5)
    xlabel('wave number [1/m]','FontSize',15);
    mk = max(kx)/5;
else             % 2D
    figure; plot(1:nk,f,'-x','linewidth',1.5)
    mk = max(kx(1:AB(1)))/2;
    hold on; plot([AB(1) AB(1)],[0 max(w(:))],'k--','linewidth',1)
    hold on; plot([AB(2) AB(2)],[0 max(w(:))],'k--','linewidth',1)
    xlabel('wave number (boundary or Brillouin zone)','FontSize',15);
end
ylabel('frequency [Hz]','FontSize',15);
box on; grid on

% mark in grey shade the potential band gaps
for i1 = 1:size(bg,1)
    hold on; patch( [kx(1) kx(end) kx(end) kx(1)], ...
                    [bg(i1,1) bg(i1,1) bg(i1,2) bg(i1,2)], ...
                    .8*[1 1 1],'EdgeColor','none')
end

% add marks for homogenized velocities
if nargin>4
    xk = linspace(0,mk,100);
    hold on; plot(xk,(cp/2/pi)*xk,'--k')
    hold on; plot(xk,(cb/2/pi)*xk.^2,'--k')
end

set(gca,'fontsize',15);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,v] = switchRows(w,v,indi,indj)
tmp = v(:,indi(1),indj);
v(:,indi(1),indj) = v(:,indi(2),indj);
v(:,indi(2),indj) = tmp;
tmp = w(indi(1),indj);
w(indi(1),indj) = w(indi(2),indj);
w(indi(2),indj) = tmp;
end
