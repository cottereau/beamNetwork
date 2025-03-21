function [kx,ky] = linearWaveNumber(maxKx,maxKy,nk)

kx = linspace(0,maxKx,nk);
kx = [ kx maxKx*ones(1,nk-2) kx(end:-1:2) ];

ky = linspace(0,maxKy,nk);
ky = [ zeros(1,nk-1) ky ky(end-1:-1:2) ];
