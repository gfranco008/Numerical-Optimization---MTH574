% MTH 574 - Numerical Optimization
% Gustavo Franco Reynoso 

function p = minsk(B, gradfk, dk)
% Minimize the sk using F(s) = gradfk'*s + 0.5*s'*B*s

% B       - approximate Hessian
% gradF   - function handle of gradient
% dk      - radius of circle constraint

global d
d = dk;
% Objective function
fun = @(s) gradfk'*s + 0.5*s'*B*s;
x0 = [1;1];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
options = optimoptions('fmincon','Display','off');
% constraints
nonlcon = @circlecon;
p = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end

% Constraint
function [c,ceq] = circlecon(x)
% This gives us the circle constraint

global d
%s(1)^2 + s(2)^2 <= dk;
c = (x(1))^2 + (x(2))^2 - d^2;
ceq = [];
end

