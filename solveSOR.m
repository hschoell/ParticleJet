function [x_1,res_end] = solveSOR(aw,ae,an,as,ap,rhs,x_1,nx,ny)
% This function solves a linear system of equations using
% Successive-Over-Relaxation. It takes the relevant diagonals of A and the
% position of which in d. Also the right hand side and an initial guess.

lim = 1e-7; % residual limit
maxIt = 100000; % iterational limit
omega = 1.5; % relaxation factor

res = [1];
% Initialize x
x_s=zeros(nx,ny);

% SOR
while res(end)>lim && length(res)<maxIt
    x_0=x_1;
    for i = 1:nx
        for j = 1:ny
            b= i+nx*(j-1);           
            x_s(i,j)=-1/ap(b)*(as(b)*x_s(i,1+mod(j-2,ny))...
                + aw(b)*x_s(1+mod(i-2,nx),j)...
                + ae(b)*x_1(1+mod(i,nx),j)...
                + an(b)*x_1(i,1+mod(j,ny)) + rhs(b));
        end
    end
    x_1 = x_1 + omega*(x_s-x_1);
    res = [res norm(x_1-x_0)];
end
res_end = res(end);

% figure(1)
% semilogy(1:length(res),res,'LineWidth',2)
% title('Residuals');
% xlabel('Iterations');
% ylabel('\epsilon_k');
end