function x_1 = solveCG(A,rhs,x_1)
% This function solves a pressure-poisson equation using the
% Conjugate-Gradient-Method. A is the operator matrix, rhs the
% right-hand-side and x_1 is the intitial guess in vector form. Also we're
% ignoring whether its spd or not

%check whether positive definite
K = full(A);
[~,p] = chol(A) % 0 if pd 
sym = issymmetric(A)
% if p 
%     A = [zeros(size(A)) A; A' zeros(size(A))];
%     rhs = [rhs;zeros(size(rhs))];
%     x_1 = [zeros(size(x_1));x_1];
% end
% 
% % Check again
% K = full(A);
% [~,p] = chol(A) % 0 if pd 

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
% x_1=flip(x_1);

% length(rho)
% x_1 = x_1 - (1/(length(x_1)))*sum(x_1);
% sum(x_1)
% figure(2)
% clf
% semilogy(1:length(rho),rho,'LineWidth',2)
% title('Residuals');
% xlabel('Iterations');
% ylabel('\epsilon_k');
end