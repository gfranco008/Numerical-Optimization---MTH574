% MTH 574 - Numerical Optimization
% Gustavo Franco Reynoso 
function hessfk = subproblem(dk, gradfk, hessfk, N)
% This  solves for the approximation of the Hessian that is ositive
% definite using:
% 
% dk      - Radius of the circle
% gradfk  - Gradient at xk
% Hessfk  - Hessian at xk
% N       - Number of iterations
% ------------------------------------------------------------------

% Determine if the hessian is pos def with cholesky decomposition
[R,p] = chol(hessfk);
lambda = eig(hessfk);
% Select smallest lambda greater than 0 or lambda = 1
lambda = min(lambda(eig(hessfk)>0));
if isempty(lambda) == 1
    lambda = 1;
end
n = 0;

if p~=0
    while p~=0 || n<N
        % Update B (approximate hessian)
        B = hessfk + lambda.*eye(size(hessfk));
        % check if it is pos def
        [R,p] = chol(B);
        %keyboard
        % If you cant decompose it then update again
        if isempty(R) == 1
            lambda = lambda*0.25;
            B = hessfk + lambda*eye(size(hessfk));
            [R,p] = chol(B);
            if isempty(R) == 1
                break;
            end
        else
            % update lambda
            plam = -B\gradfk;            
            qlam = R'\plam;
            lambda = lambda + ((norm(plam)/norm(qlam)).^2)*(norm(plam)-dk)/dk;
            n = n+1;
            if p == 0 && n >=2
                break;
            end
        end
    end
    hessfk = B;
elseif p == 0
    hessfk;
end