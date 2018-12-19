% Function to calculate del(del p) using 
function [ap,ae,aw,an,as,rhs] = eq17(uStar,vStar,dx,dy,dt,nx,ny);

% Number of unknowns
args = nx*ny;

% Initializing the RHS including the boundaries (No BCs needed, 
% since p = nx*ny but u = nx+1*ny and v
rhs = zeros(args,1);

for j = 1:ny %Iterating over y
    for i = 1:nx %Iterating over x
        rhs(i+(j-1)*nx) = 1/dt.*((uStar(i+1,j)-uStar(i,j))/dx...
            +(vStar(i,j+1)-vStar(i,j))/dy);       
    end
end

% Calculating the diagonals of the operator matrix
aw = 1/dx.^2*ones(args,1);
ae = 1/dx.^2*ones(args,1);
an = 1/dy.^2*ones(args,1);
as = 1/dy.^2*ones(args,1);
ap = -2*(1/dx.^2+1/dy.^2)*ones(args,1);
% Correct diagonal values at boundaries (zerogradient everywhere)
ap(nx+1:nx:nx*(ny-1)-nx+1)=-1/dx.^2-2/dy.^2; %western boundary
ap(2*nx:nx:nx*(ny-1))=-1/dx.^2-2/dy.^2; %eastern boundary
ap(2:nx-1)=-2/dx.^2-1/dy.^2; %south boundary 
ap(nx*ny-nx+2:nx*ny-1)=-2/dx.^2-1/dy.^2; %north boundary

ap([1,nx,nx*ny-nx+1,nx*ny])=-1/dx.^2-1/dy.^2;


for m=1:nx
    for n = 1:ny
        d = m+nx*(n-1);
        if m == 1
            aw(d) = 0;
        elseif m == nx
            ae(d) = 0;
        end
        if n == 1
            as(d) = 0;
        elseif n == ny
            an(d) = 0;
        end
    end
end

end