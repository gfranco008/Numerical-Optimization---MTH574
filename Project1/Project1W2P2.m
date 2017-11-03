%Modify the codes to find the Maximums 

% MTH 574 - Numerical Optimization
% Gustavo Franco Reynoso 
clear all; close all;

g = @(x,y) -exp(-5*(x.^2 + y.^2));
dg = @(x,y) [10*x*exp(-5*(x.^2 + y.^2));10*y*exp(-5*(x.^2 + y.^2))];
hg = @(x,y) [(100*x.^2-10)*exp(-5*(x.^2 + y.^2)),-100*x*y*exp(-5*(x.^2 + y.^2));...
    -100*x*y*exp(-5*(x.^2 + y.^2)),(100*y.^2 -10)*exp(-5*(x.^2 + y.^2))];

K = [];
Fk = [];
N = [];
Pk = [];

    [k,fk] = linesearch_steepestascentmax(g,dg)
    [n, pk] = linesearch_newtonmax(g,dg,hg)
  
    