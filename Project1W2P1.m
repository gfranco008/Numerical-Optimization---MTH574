% Test the code and see if they can handle these functions easily.
%MTH 574 - Numerical Optimization
% Gustavo Franco Reynoso 


clear all; close all;


f = @(x,y) x.^2 + y.^2;
g = @(x,y) -exp(-5*(x.^2 + y.^2));
h = @(x,y) 2*x + x.^2 - 2*(y.^2);
%This function explodes with steepest descent because you have a saddle function. When you look
%at it from one side, it can be a minimum, but when you look at it from
%another side it becomes a maximum so it creates an issue.

F = {f;g;h};

df = @(x,y) [2*x; 2*y];
dg = @(x,y) [10*x*exp(-5*(x.^2 + y.^2));10*y*exp(-5*(x.^2 + y.^2))];
dh = @(x,y) [2 + 2*x; -4*y];
%

DF = {df;dg;dh};

hf = @(x,y) [2,0;0,2];
hg = @(x,y) [(100*x.^2-10)*exp(-5*(x.^2 + y.^2)),-100*x*y*exp(-5*(x.^2 + y.^2));...
    -100*x*y*exp(-5*(x.^2 + y.^2)),(100*y.^2 -10)*exp(-5*(x.^2 + y.^2))];
%Both methods work for most of the problems. The only issue is the line
%search newton method for the second problem. When calculating pk you
%devide the hessian over the gradient and then multiply it by negative one.
%The problem arises when you try to normalize it and find the norm it gives
%you a positive 1 which is completely opposite direction of where we should
%go to find the min. The Hessian is supposed to be positive definite.
hh = @(x,y) [2,0;0,-4];
    
Hf = {hf,hg,hh};
K = [];
Fk = [];
N = [];
Pk = [];
for i = 1:length(F)
    [k,fk] = linesearch_steepestdescent(F{i}, DF{i})
    [n, pk] = linesearch_newton(F{i},DF{i},Hf{i})
    K = [K;k];
    Fk = [Fk;fk];
    N = [N;n];
    Pk = [Pk; pk];
end
    