% Create a hybrid code where steepest descent is used at a point where the
% Hessian is not positive definite and Newton is used at a point where
% Hessian is positive definite. Test your code with the Rosenbrock
% function.

% MTH 574 - Numerical Optimization
% Gustavo Franco Reynoso 

clear all; close all;

f = @(x,y) (1-x).^2 + 100*(y - x.^2).^2;
df = @(x,y) [2*(200*x.^3 - 200*x*y + x - 1); 200*(y - x.^2)];
hf = @(x,y) [1200*x.^2 - 400*y + 2, -400*x;-400*x, 200];

[k, fk] = linesearch_hybrid(f,df,hf);
    