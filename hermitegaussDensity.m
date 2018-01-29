function dens = hermitegaussDensity(a,b,phi,n,m,min,max,type)
% PLOTPROBABILITYDENSITY
% x=min:max y=min:max

% constants
B = 2/(n-0.5);

% Calculating normalization constant
normx = sqrt(1/(sqrt(pi) * 2^n * factorial(n)) * a);
normy = sqrt(1/(sqrt(pi) * 2^m * factorial(m)) * b);

% Set up Cartesian grid
[X,Y] = meshgrid(min:max,min:max);

% Calculate radial function
syms x y;
% following found hypergeom func not reliable 

hgx = poly2sym(hermite(n),x);
hgy = poly2sym(hermite(m),y);
hgx = exp(-x^2/2) * hgx;
hgy = exp(-y^2/2) * hgy;

% Discretize
hgxy = hgx * hgy;
hgfunc = matlabFunction(hgxy);
hgmat = hgfunc(a*X, b*Y);

if(strcmp(type,'pdensity'))
 
  dens = hgmat.*hgmat;
 
elseif(strcmp(type,'real'))
 
  dens = hgmat.*cos(m*(atan2(Y,X)+phi));
  prob = sum(sum(dens.*dens));
  dens = dens/sqrt(prob);
  figure, mesh(X,Y,dens);
  figure, surf(X,Y,dens);
  figure, imagesc(dens);
   
elseif(strcmp(type,'imag'))
 
  dens = hgmat.*sin(m*(atan2(Y,X)+phi));
  prob = sum(sum(dens.*dens));
  dens = dens/sqrt(prob);
  %figure, mesh(X,Y,dens);
 
end;

