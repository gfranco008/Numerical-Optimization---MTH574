% MTH 574 - Numerical Optimization
% Gustavo Franco Reynoso 
function alpha = WolfeCon(f,gradf,xk,pk)
% implementation of strong wolfe line search
% 

% f       - function handle of objective function
% gradF   - function handle of gradient
% x       - current iterate
% pk      - search direction
% alphaLo - alpha at lower objective function value for zoom
% alphaHi - alpha at higher objectvie function value for zoom

% Define initial conditions 
alpha0 = 0;
alpha1 = 0.025; %a(i-1)
alphamax = 100;
alpha = alpha1; %a(i)
c1 = 1e-4;
c2 = 0.25;
i = 1;

while 1
     % Update xknew  and define Phi, Phi'
    xc = xk + alpha*pk;
    Of = f(xc); %Phi(alpha)
    POf = gradf(xc)'*pk; %Phi'(alpha)

    if and(or(Of > f(xk) + c1*alpha*gradf(xk)'*pk, Of >= f(xk+alpha1*pk)), i>1)
      
        alpha = Zoom(alpha1,alpha,f,gradf,xk,pk,c1,c2);
        break;
        
    elseif abs(POf) <= -c2*gradf(xk)'*pk
        
        break;
        
    elseif POf >= 0
       
        alpha = Zoom(alpha,alpha1,f,gradf,xk,pk,c1,c2);
        break;
    end
    alpha1 = alpha;
    alpha = min(alphamax, alpha1 + alpha1*2*rand(1));
    i = i+1;
    xk = xc;
end
end

function alpha = Zoom(alphaLo, alphaHi,f,gradf,xk,pk,c1,c2)
while 1
    alpha = (alphaLo+alphaHi)/2;
    % Update xknew  and define Phi, Phi'
    xc    = xk + alpha*pk;
    Of = f(xc);
    POf = gradf(xc)'*pk;
    POf0 = gradf(xk)'*pk;
    
    % check if current iterate violates sufficient decrease
    if Of > f(xk) + c1*alpha*POf0 || Of >= f(xk + alphaLo*pk)
        alphaHi = alpha;
      
    else
        % current iterate has sufficient decrease?
        if abs(POf)<=-c2*POf0
           
            break;
        end
        % are we behind the minimum?
        if POf*(alphaHi - alphaLo) >= 0 
            %pause
            alphaHi = alphaLo;
            alphaLo = alpha;
        end
        
    end
    xk = xc;
   if (abs((alphaLo-alphaHi))< 10^-8)
       alpha = alphaLo;
       break;
   end
end

end