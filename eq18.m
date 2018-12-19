function [u, v] = eq18(p, uStar, vStar,dt,dx,dy)
% This function calculates dp based upon equation 18. BCs are 0 at the
% boundaries. 
nx = length(p(:,1));
ny = length(p(1,:));
dpx = zeros(nx+1,ny);
dpy = zeros(nx,ny+1);
v = zeros(size(vStar));

for i = 2:nx
    for j = 1:ny
        dpx(i,j) = (p(i,j) - p(i-1,j))/dx;
    end
end

for i = 1:nx
    for j = 2:ny
        dpy(i,j) = (p(i,j) - p(i,j-1))/dy;
    end
end

u = uStar - dt.*(dpx);
% v(2:end-1,:) = vStar(2:end-1,:) - dt.*(dpy(2:end-1,:));
v = vStar - dt.*(dpy);
end