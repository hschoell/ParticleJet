function duStarStar = eq14(H_0, H_1, u1, d1, d2, dt, Re)
% This function calculates delta u star star using eqn 14. 

% Determine Number of unknowns
nx = length(u1(:,1));
ny = length(u1(1,:));
args = nx*ny;

% Set up the right-hand side. Again; ignore the boundaries; lexiographic
% ordering
rhs = zeros(args,1);
for j = 2:ny-1 % Iterate over y
    for i = 2:nx-1 %Iterate over x
        rhs(i+(j-1)*nx) = dt/2*(3*H_1(i,j) - H_0(i,j))...
            +dt/(Re)*(1/d1^2*(u1(i+1,j)-2*u1(i,j) + u1(i-1,j))...
            +1/d2^2*(u1(i,j+1)-2*u1(i,j)+u1(i,j-1)));
    end
end

% Care for outlet boundary conditions
for j = 2:ny-1
    rhs(nx+(j-1)*nx) = dt/2*(3*H_1(end,j) - H_0(end,j))...
            +dt/Re*(1/d1^2*(-u1(i,j) + u1(i-1,j))...
            +1/d2^2*(-u1(i,j)+u1(i,j-1)));
end

% % For plotting purposes only (???)
% for k = 1:nx
%     for l = 1:ny
%         b = k+nx*(l-1);
%         rhs1(k,l) = rhs(b);
%     end
% end

% Set up the operator on the left-hand-side
aw = -dt/(2*d1^2*Re)*ones(args,1);
ap = (1+dt/(d1^2*Re))*ones(args,1);
ae = -dt/(2*d1^2*Re)*ones(args,1);

for m=1:nx
    for n = 1:ny
        d = m+nx*(n-1);
        if m == 1% &&  n ~= 1 % Don't really understand why
            aw(d) = 0;
        elseif m == nx% && n ~=ny
            ae(d) = 0;
        end
    end
end

% Care for Outlet Boundaries:
ap(2*nx:nx:end-nx) = 1+dt/(d1^2*Re);

% Solve using Thomas-Algorithm
duStarStar_vec = thomas(aw,ap,ae,rhs,args,nx);   

% Put back into matrix form
duStarStar = reshape(duStarStar_vec,[nx,ny]);
end