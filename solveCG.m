function x_1 = solveCG(A,rhs,x_1)
% This function solves a pressure-poisson equation using the
% Conjugate-Gradient-Method. A is the operator matrix, rhs the
% right-hand-side and x_1 is the intitial guess in vector form

%check whether positive definite
K = full(A);
[~,p] = chol(A) % 0 if pd 

lim = 1e-7; % residual limit
maxIt = 1000; % iterational limit

% Initialize x
d = zeros(length(x_1),1);

% Conjugate Gradient
r = rhs - A*x_1;
rho = [norm(r)^2];
while sqrt(rho(end))>lim && length(rho)<maxIt
    x_0=x_1;
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
x_1=flip(x_1);

nx = sqrt(length(x_1));
ny = sqrt(length(x_1));
figure(1)
clf
plot((nx*ny/2):1:((ny+2)*nx)/2,x_1((nx*ny/2):1:((ny+2)*nx)/2));
xlabel('x');
ylabel('\phi');
title(['Values of \phi at y=0.5 With Conjugate Gradient']);

length(rho)
x_1 = x_1 - (1/(length(x_1)))*sum(x_1);
sum(x_1)
figure(2)
clf
semilogy(1:length(rho),rho,'LineWidth',2)
title('Residuals');
xlabel('Iterations');
ylabel('\epsilon_k');
end