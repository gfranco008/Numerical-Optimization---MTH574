% Vanilla simulation of trust region based on linear model
% Alfa Heryudono MTH574/474

% Start fresh
clear all; close all;

% Define objective function
f = @(x,y) 100*(y-x.^2).^2 + (1-x).^2;

% Define gradient of objective function
gradf = @(x,y) [ -400*x.*(y-x.^2)-2*(1-x) ; 200*(y-x.^2) ];

% Set the intial starting point
xk = [0;1];

% Set Algorithm 4.1 parameters
dhat = 1; % must be greater than 0
dk = 0.6; % initial delta, between 0 and dhat
eta = 1e-2; % [0,0.25)

% Maximum Iteration
MaxIter = 15;

% Surf and Contour plot of the objection function
figure('Position',[30 100 1200 500])
N = 60;
x = linspace(-1.25,1.25,N);
[X,Y] = meshgrid(x,x);
Z = f(X,Y);

subplot(1,2,1)
surf(X,Y,Z);
xlabel('x'); ylabel('y'); zlabel('f(x,y)')
hold on;

subplot(1,2,2)
contour(X,Y,Z,N); hold on;
th = linspace(0,2*pi);unitcirc = exp(1i*th);
h = plot(NaN,NaN,'-r');
xlabel('x'); ylabel('y');

% Print iteration header
fk = f(xk(1),xk(2)); k = 0;
fprintf('iter                 xk                    fk                  deltak             rhok\n')
fprintf('-----------------------------------------------------------------------------------------\n')
fprintf('%3d        %3.4e    %3.4e       %3.4e        %3.4e         %3.4e\n',k,xk(1),xk(2),fk,dk,NaN)

% Start the iteration
KeepIterate = true;
while KeepIterate == true

    % plot xk and trust region radius
    subplot(1,2,2)
    plot(xk(1),xk(2),'o','MarkerFaceColor','r'); hold on;
    c = xk(1)+1i*xk(2) + dk*unitcirc;
    set(h,'xdata',real(c),'ydata',imag(c)); drawnow
    pause(0.5)

    % Choose the direction vector pk using steepest descent
    gradfk = gradf(xk(1),xk(2));
    pk = -dk*gradfk/norm(gradfk);

    % -------------- Algorithm 4.1 -------------------
    fk = f(xk(1),xk(2)); fkpk = f(xk(1)+pk(1),xk(2)+pk(2));
    mk0 = fk; mkpk = fk + gradfk.'*pk;
    rhok = (fk - fkpk)/(mk0 - mkpk);

    if rhok < 0.25
        dk = 0.25*dk;
    else
        if rhok > 0.75 && abs(norm(pk)-dk) < eps
            dk = min(2*dk,dhat);
        end
    end
    if rhok > eta
        % Step xk to xk + pk
        xknew = xk+pk;

        % Test whether the function value is decreasing or not.
        % One of many criteria to stop
        fktest = f(xknew(1),xknew(2));
        if fktest < fk
            % Plot the xknew
            plot(xknew(1),xknew(2),'o','MarkerFaceColor','r')

            % Plot the pk vector (arrow)
            quiver(xk(1),xk(2),pk(1),pk(2),1,'k','Linewidth',2);
            drawnow;

            % Update xk and fk
            k = k+1; xk = xknew; fk = fktest;
        else
            KeepIterate = false;
        end
    else
        k = k + 1;
    end

    if k > MaxIter
        KeepIterate = false;
    end
    % Print the iteration
    fprintf('%3d        %3.4e    %3.4e       %3.4e        %3.4e         %3.4e\n',k,xk(1),xk(2),fk,dk,rhok)
end
