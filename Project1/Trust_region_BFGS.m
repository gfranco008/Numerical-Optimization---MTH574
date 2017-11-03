% Trust region BFGS
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
    
    
    % Set Algorithm 4.1 parameters
    dhat = 10; % must be greater than 0
    dk = 0.6; % initial delta, between 0 and dhat
    eta = 1e-2; % [0,0.25)
    
    % Set Algorithm 4.3 parameters
    n = 200;
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
    fprintf('                                  \n')
    fprintf('iter                 xk                    fk                  deltak             rhok\n')
    fprintf('-----------------------------------------------------------------------------------------\n')
    fprintf('%3d        %3.4e    %3.4e       %3.4e        %3.4e         %3.4e\n',k,xk(1),xk(2),fk,dk,NaN)
    
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
        
        % Choose the direction vector pk using steepest descent
        
       
        % -------------- Algorithm 6.1 -------------------
       
       gradfk = gradf(xk(1),xk(2));
        pk = -(B)*gradfk;
        pk = dk*pk/norm(pk);
        alpha = WolfeConBFGS(f,gradf,xk,pk); 
        % -------------- Algorithm 4.1 -------------------
        fk = f(xk(1),xk(2)); fkpk = f(xk(1)+pk(1),xk(2)+pk(2));
        mk0 = fk; mkpk = fk + gradfk.'*pk ;
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
            xknew = xk+alpha*pk;
            
            % Test whether the function value is decreasing or not.
            % One of many criteria to stop
            fktest = f(xknew(1),xknew(2));
            if fktest < fk
                % Plot the xknew
                plot(xknew(1),xknew(2),'o','MarkerFaceColor','r')
                
                % Plot the pk vector (arrow)
                quiver(xk(1),xk(2),pk(1),pk(2),1,'k','Linewidth',2);
                drawnow;
                
                % Update xk and fk and algorithm 6.1
                sk = xknew - xk;
                yk = gradf(xknew(1),xknew(2))-gradf(xk(1),xk(2));
                Bnew = (eye(size(B)) - rhok*sk*yk')*inv(B)*(eye(size(B)) - rhok*yk*sk') + rhok*sk*sk';
                k = k+1; xk = xknew; fk = fktest; B = (Bnew); 
            else
                j = j+1;
                if j > 15
                    k = k+1;
                end
            end
        else
            k = k+1;
        end
        
        if k > MaxIter
            break;
        end
        
        % Print the iteration
        fprintf('%3d        %3.4e    %3.4e       %3.4e        %3.4e         %3.4e\n',k,xk(1),xk(2),fk,dk,rhok)
 
    end
end
