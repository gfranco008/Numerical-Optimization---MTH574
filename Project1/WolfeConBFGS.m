% MTH 574 - Numerical Optimization
% Gustavo Franco Reynoso 
function alpha = WolfeConFBGS(f,gradf,xk,pk)
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
alpha1 = 1; %a(i-1)
alphamax = 100;
alpha = alpha1; %a(i)
c1 = 1e-4;
c2 = 0.9;
i = 1;

while 1
    % Update xknew  and define Phi, Phi'
    xc = xk + alpha*pk;
    Of = f(xc(1),xc(2)); %Phi(alpha)
    POf = gradf(xc(1),xc(2))'*pk; %Phi'(alpha)

    if and(or(Of > f(xk(1),xk(2)) + c1*alpha*gradf(xk(1),xk(2))'*pk, Of >= f(xk(1)+alpha1*pk(1),xk(2)+alpha1*pk(2))), i>1)
     
        alpha = Zoom(alpha1,alpha,f,gradf,xk,pk,c1,c2);
        break;
        
    elseif abs(POf) <= -c2*gradf(xk(1),xk(2))'*pk
        
        
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
    Of = f(xc(1),xc(2));
    POf = gradf(xc(1),xc(2))'*pk;
    POf0 = gradf(xk(1),xk(2))'*pk;
    if Of > f(xk(1),xk(2)) + c1*alpha*POf0 || Of >= f(xk(1) + alphaLo*pk(1), xk(2) + alphaLo*pk(2))
        alphaHi = alpha;
       
    else
        if abs(POf)<=-c2*POf0
            
            break;
        end
        if POf*(alphaHi - alphaLo) >= 0
          
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