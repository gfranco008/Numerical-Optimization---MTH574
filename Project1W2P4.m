% Vanilla simulation of a line search algorithm based on steepest descent.
% Alfa Heryudono MTH574/474

% Start fresh
clear all; close all;

% Define objective function
f = @(x,y) 100*(y-x.^2).^2 + (1-x).^2;

% Define gradient of objective function
gradf = @(x,y) [ -400*x.*(y-x.^2)-2*(1-x) ; 200*(y-x.^2) ];

% Set the intial starting point
xk = [0;1];

% Set the step-length alpha (fixed for now)
alpha = 0.0025;
delta = 0.5;
gamma = 1e-4;
% Surf and contour plot of the objection function
figure('Position',[30 100 1200 500])
x = linspace(-1.25,1.25,60); y = linspace(-0.4,1.2,60);
[X,Y]=meshgrid(x,y);
Z = f(X,Y);
subplot(1,2,1)
surf(X,Y,Z)
subplot(1,2,2)
contour(X,Y,Z,60); hold on;
plot(1,1,'o','MarkerFaceColor','r','MarkerSize',4)
hold on

% Plot initial point
plot(xk(1),xk(2),'o','MarkerFaceColor','r')

% Print iteration header
fk = f(xk(1),xk(2)); k = 0;
fprintf('iter                   xk                          fk \n')
fprintf('------------------------------------------------------------------\n')
fprintf('%3d        %3.8e    %3.8e       %3.8e   %3.3d\n',k,xk(1),xk(2),fk,alpha)

% Start the iteration
KeepIterate = true;
while KeepIterate == true

    % Choose the direction vector pk using steepest descent
    gradfk = gradf(xk(1),xk(2));
    pk = -gradfk/norm(gradfk); %normalize the vector
    % Step xk to xk + alpha*pk
    xknew = xk+alpha*pk;
    
    % Set the step-length alpha (variable) 
         
    if (f(xknew(1),xknew(2))<=f(xk(1),xk(2))+gamma*alpha*gradf(xk(1),xk(2))'*pk)
        alpha = alpha+.0005;
        if (f(xknew(1),xknew(2))<=f(xk(1),xk(2))+gamma*alpha*gradf(xk(1),xk(2))'*pk)
        else
            alpha = alpha - .005;
        end
    else
        alpha = alpha*delta;
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
        if and(or(alpha > 0.05, alpha < 0.00009),k>20)
            alpha = 0.025;
        end
        fprintf('%3d        %3.8e    %3.8e       %3.8e   %3.3d\n',k,xk(1),xk(2),fk,alpha)
 else
        KeepIterate = false;
    end

end
