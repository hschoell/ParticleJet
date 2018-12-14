% Create uniform mesh
clear all; close all; clc;

% Parameters
dt = 1e-2;
maxt = 1;
t = 0;
Re = 100;
St = Re*1e-2;
R = 0.05;

nx=128;                        % Number of grid points in x
ny=64;                         % Number of grid points in y
x=linspace(0,2,nx+1);          % x-Grid (location of cell faces)
y=linspace(0,1,ny+1);          % y-Grid (location of cell faces)
xm=x(1:end-1)+(x(2)-x(1))/2;   % x-Grid (location of cell centers)
ym=y(1:end-1)+(y(2)-y(1))/2;   % y-Grid (location of cell centers)
dx=xm(2)-xm(1);                % Grid spacing in x
dy=ym(2)-ym(1);                % Grid spacing in y

%setting up initial conditions, zero everywhere
u = zeros(nx+1,ny); %store velocities at face centers
v = zeros(nx,ny+1); %store velocities at face centers
p = zeros(nx,ny);
dp = zeros(nx,ny);

u(1,:) = exp(-(ym - 0.5).^2./(R^2)); %initialize u with jet conditions

% Convective terms for u and v
Hu_0 = zeros(nx+1,ny);
Hu_1 = zeros(nx+1,ny);
Hv_0 = zeros(nx,ny+1);
Hv_1 = zeros(nx,ny+1);

while t<maxt
    % Step 1
    % setting up H
    Hu_1 = convec(u,v,dx,dy);
    Hv_1 = convec(v',u',dy,dx)';    
    
    % Solving eqn 14 for duStarStar
    duStarStar = eq14(Hu_0, Hu_1, u, dx, dy, dt, Re);
    dvStarStar = eq14(Hv_0, Hv_1, v, dx, dy, dt, Re);
    
    % Solving eqn 15 for duStar
    duStar = eq15(duStarStar, dy, dt, Re);
    dvStar = eq15(dvStarStar, dy, dt, Re);
    
    uStar = u + duStar;
    vStar = v + dvStar;
    
    %Solving eqn 17 for del*(del p)
    [ap,ae,aw,an,as,rhs] = eq17(uStar,vStar,u,v,dx,dy,dt,nx,ny);
       
    p_0 = solveSOR(aw,ae,an,as,ap,rhs,p,nx,ny);
    
    for i = 1:nx
        for j = 1:ny
            b = i+nx*(j-1);
            p(i,j) = p_0(b);
        end
    end
        
    for i = 2:nx
        for j = 2:ny
            dp(i,j) = (p(i,j) - p(i-1,j))/dx + (p(i,j) - p(i,j-1))/dy;
        end
    end
    
    u = uStar - dt*(dp(:,1:ny));
    v = vStar - dt*(dp(1:nx,:));
    
    Hu_0 = Hu_1;
    Hv_0 = Hv_1;
    t = t + dt;
    
    surf(u)
end


