function dens = laguerreDensity(a,phi,n,m,min,max,type)
% x=min:max y=min:max

% constants
B = 2/(n-0.5);

% Set up Cartesian grid
[X,Y] = meshgrid(min:max,min:max);

% Calculate radial function
syms x r;
%laguerre = hypergeom(-n+abs(m)+1,2*abs(m)+1,x);
%laguerre = subs(laguerre,x,B*r);
%radial = (B*r)^abs(m) * exp(-B*r/2) * laguerre
laguerre = poly2sym(LaguerreGen(n-abs(m)-1,2*abs(m)),r);
laguerre = exp(-r/2) * r^m * laguerre;

% Discretize
%RADIAL = subs(radial,r,a*sqrt(X.^2+Y.^2)); % put in a * sqrt(X.^2+Y.^2) for scale factor
lagfunc = matlabFunction(laguerre);
lagmat = lagfunc(a*sqrt(X.^2+Y.^2));

if(strcmp(type,'pdensity'))
 
  dens = lagmat.*lagmat;
 
elseif(strcmp(type,'real'))
 
  dens = lagmat.*cos(m*(atan2(Y,X)+phi));
  prob = sum(sum(dens.*dens));
  dens = dens/sqrt(prob);
%  figure, mesh(X,Y,dens);
%  figure, surf(X,Y,dens);
%  figure, imshow(dens);
 
elseif(strcmp(type,'imag'))
 
  dens = lagmat.*sin(m*(atan2(Y,X)+phi));
  prob = sum(sum(dens.*dens));
  dens = dens/sqrt(prob);
  %figure, mesh(X,Y,dens);
 
end;


