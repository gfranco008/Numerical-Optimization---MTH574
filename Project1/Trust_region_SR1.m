% Trust region SR1
% MTH 574 - Numerical Optimization
% Gustavo Franco Reynoso 

% Start fresh
clear all; close all; clc

% Define objective function
f = @(x,y) 100*(y-x.^2).^2 + (1-x).^2;
j = @(x,y) x.^2 + y.^2;
g = @(x,y) -exp(-5*(x.^2 + y.^2)); % This is a saddle and presents errors
h = @(x,y) 2*x + x.^2 - 2*(y.^2);

F = {f;g;h;j};

% Define gradient of objective function
df = @(x,y) [ -400*x.*(y-x.^2)-2*(1-x) ; 200*(y-x.^2) ];
dj = @(x,y) [2*x; 2*y];
dg = @(x,y) [10*x*exp(-5*(x.^2 + y.^2));10*y*exp(-5*(x.^2 + y.^2))];
%The only issue is the third problem. When calculating pk you
%devide the hessian over the gradient and then multiply it by negative one.
%The problem arises when you try to normalize it and find the norm it gives
%you a positive 1 which is completely opposite direction of where we should
%go to find the min. The Hessian is supposed to be positive definite.
dh = @(x,y) [2 + 2*x; -4*y];

gradF = {df;dg;dh;dj};


for i = 1:length(F)
    f = F{i};
    gradf = gradF{i};
    
    % Set the intial starting point
    xk = [0;1];
    xknew = xk +.025;
    
    % Set Algorithm 6.2 parameters
    dk = 0.6; % initial delta, between 0 and dhat
    nu = 10^-4; %(0,10^-3)
    r = 0.5; %(0,1)
    
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
    fprintf('iter                 xk                    fk                  deltak         \n')
    fprintf('-----------------------------------------------------------------------------------------\n')
    fprintf('%3d        %3.4e    %3.4e       %3.4e        %3.4e    \n',k,xk(1),xk(2),fk,dk)
    
    % Start the iteration
    B =eye(length(xk));
    j=0;
    
    KeepIterate = true;
    while norm(gradf(xk(1),xk(2))) > 10^-3 %KeepIterate == true
        
        % plot xk and trust region radius
        subplot(1,2,2)
        plot(xk(1),xk(2),'o','MarkerFaceColor','r'); hold on;
        c = xk(1)+1i*xk(2) + dk*unitcirc;
        set(h,'xdata',real(c),'ydata',imag(c)); drawnow
        pause(0.5)
        
        % Evaluate f and df
        fk = f(xk(1),xk(2));
        gradfk = gradf(xk(1),xk(2));
        
        
        % -------------- Algorithm 6.2 -------------------
        %solve sk
        
        %ss  = @(s) gradfk'*s + 0.5*s'*B*s ;% subject to the circle region
        % p1^2 + p2^2 <= delta
        sk = minsk(B, gradfk, dk);
        % Step xk to xk + pk
        xs = xk+sk;
        yk = gradf(xs(1),xs(2))-gradfk;
        ared = fk - f(xs(1),xs(2));
        pred = -(gradfk'*sk + 0.5*sk'*B*sk);
        if ared/pred >nu
            xknew = xk+sk;
        else
            xknew = xk;
        end
        if ared/pred >0.75
            if norm(sk) <= 0.8*dk
            else
                dk = 2*dk;
            end
        elseif (0.1 <= ared/pred )&& (ared/pred <=0.75)
        else
            dk = 0.5*dk;
        end
        if abs(sk'*(yk - B*sk))>= r*norm(sk)*norm((yk - B*sk))
            B = B + ((sk - B*yk)*(sk - B*yk)')/((sk - B*yk)'*yk);
        else
            Bnew = B;
        end
        
        % ------------------------------------------------
        % Test whether the function value is decreasing or not.
        % One of many criteria to stop
        fktest = f(xknew(1),xknew(2));
        if fktest < fk
            % Plot the xknew
            plot(xknew(1),xknew(2),'o','MarkerFaceColor','r')
            
            % Plot the pk vector (arrow)
            quiver(xk(1),xk(2),sk(1),sk(2),1,'k','Linewidth',2);
            drawnow;
            
            % Update xk and fk
            k = k+1; xk = xknew; fk = fktest; %B = (Bnew);
        else
            j = j+1;
            if j > 15
                k = k+1;
            end
        end
        
        if k > MaxIter
            break;
        end
        
        % Print the iteration
        fprintf('%3d        %3.4e    %3.4e       %3.4e        %3.4e       \n',k,xk(1),xk(2),fk,dk)
        
    end
end

