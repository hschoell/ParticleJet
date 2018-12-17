% Function to calculate del(del p) using 
function [ap,ae,aw,an,as,rhs] = eq17(uStar,vStar,dx,dy,dt,nx,ny);

% Number of unknowns
args = nx*ny;

% Initializing the RHS including the boundaries (No BCs needed, since p =
% nx*ny but u = nx+1*ny and v
rhs = zeros(args,1);

for j = 1:ny %Iterating over y
    for i = 1:nx %Iterating over x
        rhs(i+(j-1)*nx) = 1/dt.*((uStar(i+1,j)-uStar(i,j))/dx...
            +(vStar(i,j+1)-vStar(i,j))/dy);       
    end
end

% For plotting puposes only
%     for k = 1:nx
%         for l = 1:ny
%             b = k+nx*(l-1);
%             rhs1(k,l) = rhs(b);
%         end
%     end

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
ap(nx*ny-nx+2:nx*ny-1)=-2/dx.^2-1/dy.^2; 

ap([1,nx,nx*ny-nx+1,nx*ny])=-1/dx.^2-1/dy.^2;


ap1 = reshape(ap,[nx,ny]);


for m=1:nx
    for n = 1:ny
        d = m+nx*(n-1);
        if m == 1% &&  n ~= 1 % Don't really understand why
            aw(d) = 0;
        elseif m == nx% && n ~=ny
            ae(d) = 0;
        end
        if n == 1% && m ~= nx
            as(d) = 0;
        elseif n == ny% && m ~= nx
            an(d) = 0;
        end
    end
end

% For debugging purposes???
for s=1:nx
    for t = 1:ny
        e = s+nx*(t-1);
        aw1(s,t) = aw(e);
        as1(s,t) = as(e);
        ae1(s,t) = ae(e);
        an1(s,t) = an(e);
    end
end

end