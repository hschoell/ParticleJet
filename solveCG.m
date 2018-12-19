function x_1 = solveCG(A,rhs,x_1)
% This function solves a pressure-poisson equation using the
% Conjugate-Gradient-Method. A is the operator matrix, rhs the
% right-hand-side and x_1 is the intitial guess in vector form. Also we're
% ignoring whether its spd or not

%check whether positive definite
[~,p] = chol(A); % 0 if pd 
%sym = issymmetric(A); %this was for checking purposes only
if p 
    A = -1.*A;
    rhs = -1.*rhs;
end

lim = 1e-7; % residual limit
maxIt = 1000; % iterational limit

% Initialize x
d = zeros(length(x_1),1);

% Conjugate Gradient
r = rhs - A*x_1;
rho = [norm(r)^2];
while sqrt(rho(end))>lim && length(rho)<maxIt
    if length(rho)==1
        d = r;
    else
        beta = rho(end)/rho(end-1);
        d = r + beta*d;
    end
    e = A*d;
    a = rho(end)/(d'*e);
    x_1 = x_1 + a*d;
    r = r - a*e;
    rho = [rho r'*r];
end

end