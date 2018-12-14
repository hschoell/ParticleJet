function x_1 = solveSOR(aw,ae,an,as,ap,rhs,x_1,nx,ny)
% This function solves a linear system of equations using
% Successive-Over-Relaxation. It takes the relevant diagonals of A and the
% position of which in d. Also the right hand side and an initial guess.

lim = 1e-7; % residual limit
maxIt = 10000; % iterational limit
omega = 1.5085; % relaxation factor

res = [1];
% Initialize x
x_s=zeros(nx,ny);

% SOR
while res(end)>lim && length(res)<maxIt
    x_0=x_1;
    for i = 1:nx
        for j = 1:ny
            x_s(i,j)=-1/ap(i,j)*(as(i,j)*x_s(i,1+mod(j-2,ny))...
                + aw(i,j)*x_s(1+mod(i-2,nx),j)...
                + ae(i,j)*x_1(1+mod(i,nx),j)...
                + an(i,j)*x_1(i,1+mod(j+1,ny)) + rhs(i));
        end
    end
    x_1 = x_1 + omega*(x_s-x_1);
    res = [res norm(x_1-x_0)];
end
figure(1)
semilogy(1:length(res),res,'LineWidth',2)
title('Residuals');
xlabel('Iterations');
ylabel('\epsilon_k');
end