function dens = hermitegaussDensity1D(a,n,min,max)
% PLOTPROBABILITYDENSITY
% x=min:max y=min:max

% constants
B = 2/(n-0.5);

% Calculating normalization constant
normx = sqrt(1/(sqrt(pi) * 2^n * factorial(n)) * a);

% Set up Cartesian grid
[X] = linspace(min,max,100);

% Calculate radial function
syms x r;

her = hermite(n); 
hgx = poly2sym(her,x);
hgx = exp(-x^2/2) * hgx;
figure, ezplot(hgx);



