% Line search algorithm with Wolfe conditions
% MTH 574 - Numerical Optimization
% Gustavo Franco Reynoso 

% Start fresh
clear all; close all;

% Objective function
f = @(x) x(1).^2 + x(2).^2;
g = @(x) -exp(-5*(x(1).^2 + x(2).^2));
h = @(x) sin(pi*x(1))*sin(pi*x(2));
%This function explodes with steepest descent because you have a saddle function. When you look
%at it from one side, it can be a minimum, but when you look at it from
%another side it becomes a maximum so it creates an issue.

F = {f;g;h};

% Gradient of Objective function
df = @(x) [2*x(1); 2*x(2)];
dg = @(x) [10*x(1)*exp(-5*(x(1).^2 + x(2).^2));10*x(2)*exp(-5*(x(1).^2 + x(2).^2))];
dh = @(x) [pi.*cos(pi.*x(1)).*sin(pi.*x(2));pi.*cos(pi.*x(2)).*sin(pi.*x(1))];

Df = {df;dg;dh};


fp = @(x,y) x.^2 + y.^2;
gp = @(x,y) -exp(-5*(x.^2 + y.^2));
hp = @(x,y) sin(pi*x).*sin(pi*y);

Fp ={fp,gp,hp};


for j = 1:3
    f=F{j};
    gradf = Df{j};
    
    % Surf and Contour plot of the objection function
    figure('Position',[30 100 1200 500])
    N = 60;
    x = linspace(-1,1,N);
    [X,Y] = meshgrid(x,x);
    
    fplot =Fp{j};
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
    title('Wolfe Condition');
    
    % Set initial point
    xk = [-0.25;-0.5];
    xc = [0;0];
    fk = f(xk);
    e1 = 10^-4;
    e2 = 10^-4;
    k = 0;
    fprintf(' \n')
    fprintf('                         Line Search Wolfe Conditions \n')
    fprintf('iter                   xk                          fk                      alphak\n')
    fprintf('-----------------------------------------------------------------------------------------\n')
    fprintf('%3d        %3.8e    %3.8e       %3.8e        %3.8e\n',0,xk(1),xk(2),fk, 1)
    
    % Perform interation
    while ~and(norm(xk-xc)<e1, norm(gradf(xk)-gradf(xc))<e2)
        subplot(1,3,2)
        plot(xk(1),xk(2),'o','MarkerFaceColor','r')
        
        % Choose pk
        gradfk = gradf(xk); % or you can normalize this
        pk = -gradfk;
        
        % Plot Wolfe Curve
        % 
        c = 0.25; rho = 0.95;
        alpha = linspace(0,1,121);
        xx = bsxfun(@plus,xk,bsxfun(@times,alpha,pk));
        set(h1,'Xdata',alpha,'YData',fplot(xx(1,:),xx(2,:)));
        set(h2,'Xdata',alpha,'YData',fplot(xk(1),xk(2)) + (c*gradfk.'*pk)*alpha);
        
        % Choose alphak using Wolfe's Condition
        alpha = WolfeCon(f,gradf,xk,pk) ;%
        
        set(h3,'XData',alpha,'YData',f(xk + alpha*pk));
        set(h4,'XData',alpha,'YData',f(xk) + (c*gradfk.'*pk)*alpha);
        drawnow;
        pause(0.8)
     
        subplot(1,3,2)
        quiver(xk(1),xk(2),alpha*pk(1),alpha*pk(2),1,'k','Linewidth',2);
        
        % Step xk to xk + alpha*pk
        xc = xk;
        xk = xk+alpha*pk;
        fk = f(xk); % Update function value
        k = k+1;
        
        fprintf('%3d        %3.8e    %3.8e       %3.8e       %3.4d\n',k,xk(1),xk(2),fk,alpha)
        
    end
end
