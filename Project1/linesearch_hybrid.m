% Hybrid code where steepest descent is used at a point where the Hessian
% is not positive definite and Newton is used at a point where Hessian is
% positive definite.

% MTH 574 - Numerical Optimization
% Gustavo Franco Reynoso 

% Start fresh
%close all;
function [k, fk] = linesearch_hybrid(F,df,Hf)
% Define objective function
f = F;

% Define gradient of objective function
gradf = df;

% Define hessian matrix
hessf = Hf;

% Set the intial starting point
xk = [0;0.5];

% Set the step-length alpha (fixed for now)
alpha = 0.025;

% Surf and contour plot of the objection function
figure('Position',[30 100 1200 500])
x = linspace(-1.25,1.25,60); y = linspace(-0.4,1.2,60);
[X,Y]=meshgrid(x,y);
Z = f(X,Y);
subplot(1,2,1)
surf(X,Y,Z)
subplot(1,2,2)
contour(X,Y,Z,60); hold on;
plot(1,1,'o','MarkerFaceColor','k','MarkerSize',6)
title('Line Search Hybrid');
hold on

% Plot initial point
plot(xk(1),xk(2),'o','MarkerFaceColor','r')

% Print iteration header
fk = f(xk(1),xk(2)); k = 0;
fprintf(' \n')
fprintf('                         Line Search Hybrid\n')
fprintf('iter                   xk                          fk \n')
fprintf('------------------------------------------------------------------\n')
fprintf('%3d        %3.8e    %3.8e       %3.8e\n',k,xk(1),xk(2),fk)

% Start the iteration
KeepIterate = true;
while KeepIterate == true
    
    % Choose the direction vector pk using steepest descent Newton
    gradfk = gradf(xk(1),xk(2));
    %keyboard
    hessfk = hessf(xk(1),xk(2));
    [R,p] = chol(hessfk);
    if p == 0
        pk = -(hessfk\gradfk);
        pk = pk/norm(pk);
    else
        % Choose the direction vector pk using steepest descent
        pk = -gradfk/norm(gradfk); %normalize the vector
    end
    % Step xk to xk + alpha*pk
    xknew = xk+alpha*pk;
    
    % Test whether the function value is decreasing or not.
    % One of many criteria to stop
    fktest = f(xknew(1),xknew(2));
    if fktest < fk
        % Plot the xknew
        plot(xknew(1),xknew(2),'o','MarkerFaceColor','r')
        
        % Plot the pk vector (arrow)
        quiver(xk(1),xk(2),alpha*pk(1),alpha*pk(2),1,'k','Linewidth',2);
        drawnow;
        
        % Update xk and fk and print the iteration
        k = k+1; xk = xknew; fk = fktest;
        fprintf('%3d        %3.8e    %3.8e       %3.8e\n',k,xk(1),xk(2),fk)
    else
        KeepIterate = false;
    end

end