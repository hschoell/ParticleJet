%% This is a Navier-Stokes-Simulation of a particle laden jet.
% Credit: Daniel Mayers, Christopher Tossas & Henry Schoeller (Team CREAM)
% Here, the domain is generated and the ICs are set.

clear all; close all; clc;

% Parameters
dt = 2e-2;
maxt = 25;
t = 0;
Re = 100;
St = Re*1e-2;
R = 0.05;

% Create uniform mesh
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

u(1,2:ny-1) = exp(-(ym(1,2:ny-1) - 0.5).^2./(R^2));    %initialize u with jet conditions
u_in = sum(u(1,:));                     %total incoming flow

% Convective terms for u and v
Hu_0 = zeros(nx+1,ny);
Hu_1 = zeros(nx+1,ny);
Hv_0 = zeros(nx,ny+1);
Hv_1 = zeros(nx,ny+1);

%% Advancing in time

while t<maxt
    % Step 1
    % setting up H
    Hu_1 = convec(u,v,dx,dy);
    Hv_1 = convec(v',u',dy,dx)';    
    
    % Outlet Boundaries for the H terms
    for j = 2:ny-1    % for u
        Hu_1(end,j) = -1/(4*dx)*(3*u(end,j)^2 - 2*u(end,j)*u(end-1,j) - u(end-1,j)^2);
        Hu_1(end,j) = Hu_1(end,j) - 1/(4*dy)*((u(end,j+1)+u(end,j))*(v(end,j+1)+v(end-1,j+1)) - ...
            (u(end,j)+u(end,j-1))*(v(end,j)+v(end-1,j)));
    end
    for j = 2:ny      % for v
        Hv_1(end,j) = -1/(4*dy)*(v(end,j+1)^2 + 2*v(end,j+1)*v(end,j) - ...
            2*v(end,j)*v(end,j-1) - v(end,j-1)^2);
        Hv_1(end,j) = Hv_1(end,j) - 1/(4*dx)*(2*v(end,j)*(u(end,j)+u(end,j-1)) - ...
            (v(end,j)+v(end-1,j))*(u(end,j)+u(end,j-1)));
    end
    
    % Solving eqn 14 for duStarStar
    duStarStar = eq14(Hu_0, Hu_1, u, dx, dy, dt, Re);
    dvStarStar = eq14(Hv_0, Hv_1, v, dx, dy, dt, Re);
       
    
    % Solving eqn 15 for duStar
    duStar = eq15(duStarStar, dy, dt, Re);
    dvStar = eq15(dvStarStar, dy, dt, Re);
    
    
    uStar = u;
    vStar = v;
    uStar(2:nx+1,2:ny-1) = uStar(2:nx+1,2:ny-1) + duStar(2:nx+1,2:ny-1);
    vStar(2:nx,2:ny) = vStar(2:nx,2:ny) + dvStar(2:nx,2:ny);
    
    % Make sure, mass is conserved (additively):
    u_diff = u_in - sum(uStar(nx+1,:));
    uStar(nx+1,:) = uStar(nx+1,:)+(u_diff/ny);
    sum(uStar(1,:)-uStar(nx+1,:))
    
    %Solving eqn 17 for del*(del p)
    [ap,ae,aw,an,as,rhs] = eq17(uStar,vStar,dx,dy,dt,nx,ny);
       
    %[p_0,res] = solveSOR(aw,ae,an,as,ap,rhs,p,nx,ny);
    
    B=[an ae ap aw as];
    d=[-nx -1 0 1 nx];
    A=(spdiags(B, d, nx*ny, nx*ny))';
    p_test = A\rhs;
    p_0 = reshape(p_test,[nx,ny]);
    [u, v] = eq18(p_0, uStar, vStar,dt,dx,dy);
    
    % Output the divergence of the velocity field to check solenoidality 
    max_div = max(max((u(2:end,:)-u(1:end-1,:))/dx + ...
        (v(:,2:end)-v(:,1:end-1))/dy));
    
    Hu_0 = Hu_1;
    Hv_0 = Hv_1;
    t = t + dt
    p = p_0;
    figure(1)
    contourf(u)
    figure(3)
    contourf(v)
end