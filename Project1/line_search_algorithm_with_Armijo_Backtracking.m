% Vanilla simulation of line search algorithm with Armijo Backtracking
% Alfa Heryudono MTH574/474 - Fall 2014

% Start fresh
clear all; close all;

% Objective function
f = @(x) sin(pi*x(1)).*sin(pi*x(2));

% Gradient of Objective function
gradf = @(x)[pi.*cos(pi.*x(1)).*sin(pi.*x(2));pi.*cos(pi.*x(2)).*sin(pi.*x(1))];

% Surf and Contour plot of the objection function
figure('Position',[30 100 1200 500])
N = 60;
x = linspace(-1,1,N);
[X,Y] = meshgrid(x,x);
fplot = @(x,y) sin(pi*x).*sin(pi*y);
Z = fplot(X,Y);
subplot(1,3,1)
surf(X,Y,Z);
xlabel('x'); ylabel('y'); zlabel('f(x,y)')
hold on;
subplot(1,3,2)
contour(X,Y,Z,60)
hold on
xlabel('x'); ylabel('y');
subplot(1,3,3)
h1 = plot(NaN,NaN,'-b'); hold on; h2 = plot(NaN,NaN,'--r');
h3 = plot(NaN,NaN,'*k','MarkerSize',10); h4 = plot(NaN,NaN,'*r','MarkerSize',10);
xlabel('alpha'); ylabel('f(x,y)');
ylim([min(Z(:))-0.2 max(Z(:))+0.2]);
title('Armijo Condition');

% Set initial point
xk = [-0.25;-0.5];
fk = f(xk);

fprintf('iter                   xk                          fk                      alphak\n')
fprintf('-----------------------------------------------------------------------------------------\n')
fprintf('%3d        %3.8e    %3.8e       %3.8e        %3.8e\n',0,xk(1),xk(2),fk, 1)

% Perform interation
for k=1:10
    subplot(1,3,2)
    plot(xk(1),xk(2),'o','MarkerFaceColor','r')

    % Choose pk
    gradfk = gradf(xk); % or you can normalize this
    pk = -gradfk;

    % Plot Armijo Curve
    % Armijo Constant
    c = 0.25; rho = 0.95;
    alpha = linspace(0,1,121);
    xx = bsxfun(@plus,xk,bsxfun(@times,alpha,pk));
    set(h1,'Xdata',alpha,'YData',fplot(xx(1,:),xx(2,:)))
    pause(0.5)
    set(h2,'Xdata',alpha,'YData',fplot(xk(1),xk(2)) + (c*gradfk.'*pk)*alpha);

    % Choose alphak using Armijo's backtracking
    alpha = 1; % Starting alpha
    while f(xk + alpha*pk) > f(xk) + c*alpha*gradfk.'*pk
        alpha = rho*alpha;
        set(h3,'XData',alpha,'YData',f(xk + alpha*pk));
        set(h4,'XData',alpha,'YData',f(xk) + (c*gradfk.'*pk)*alpha);
        drawnow;
        pause(0.08)
    end

    subplot(1,3,2)
    quiver(xk(1),xk(2),alpha*pk(1),alpha*pk(2),1,'k','Linewidth',2);

    % Step xk to xk + alpha*pk
    xk = xk+alpha*pk;
    fk = f(xk); % Update function value

    fprintf('%3d        %3.8e    %3.8e       %3.8e       %3.8e\n',k,xk(1),xk(2),fk,alpha)

end
