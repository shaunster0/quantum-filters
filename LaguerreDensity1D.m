function dens = LaguerreDensity1D(a,m,n,min,max)
% PLOTPROBABILITYDENSITY
% x=min:max y=min:max


% Calculate radial function
syms r;

laguerre = poly2sym(LaguerreGen(n-abs(m)-1,2*abs(m)),r);
laguerre = exp(-r/2) * r^m * laguerre;
figure, ezplot(laguerre, [0,50]);