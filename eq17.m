% Function to calculate del(del p) using 
function [A_ddp, rhs] = eq17(uStar,vStar,u,v,dx,dy,dt,nx,ny);

% Number of unknowns
args = length(u(1,:))*length(u(:,1));

% Initializing the RHS ignoring the boundaries. 
rhs = zeros(nx*ny,1);

for j = 2:length(u(1,:))-1 %Iterating over y
    for i = 2:length(v(:,1))-1 %Iterating over x
        rhs(i+(j-1)*length(v(:,1))) = 1/dt.*((uStar(i+1,j)-uStar(i,j))/dx...
            +(vStar(i,j+1)-vStar(i,j))/dy);       
    end
end


    for k = 1:128
        for l = 1:64
            b = k+128*(l-1);
            rhs1(k,l) = rhs(b);
        end
    end

aw = 1/dx.^2*ones(args,1);
ae = 1/dx.^2*ones(args,1);
an = 1/dx.^2*ones(args,1);
as = 1/dx.^2*ones(args,1);
ap = -4/dx.^2*ones(args,1);

A = diag(ap)+diag(an(1+nx/2:end-nx/2),nx)+diag(as(1+nx/2:end-nx/2),-nx) ...
    +diag(ae(1:end-1),1)+diag(aw(1:end-1),-1);
A_ddp = sparse(A);

end







