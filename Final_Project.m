%% This is a Navier-Stokes-Simulation of a particle laden jet.
% Credit: Daniel Mayers, Christopher Tossas & Henry Schoeller (Team CREAM)
% Here, the domain is generated and the ICs are set.

clear all; close all; clc;

% Parameters
CFL = 0.45;
maxt = 2;
t = 0;
Re = 100;
St = 1;
R = 0.05;
i = 1;

% Create uniform mesh
nx=128;                        % Number of grid points in x
ny=64;                         % Number of grid points in y
x=linspace(0,2,nx+1);          % x-Grid (location of cell faces)
y=linspace(0,1,ny+1);          % y-Grid (location of cell faces)
xm=x(1:end-1)+(x(2)-x(1))/2;   % x-Grid (location of cell centers)
ym=y(1:end-1)+(y(2)-y(1))/2;   % y-Grid (location of cell centers)
dx=xm(2)-xm(1);                % Grid spacing in x
dy=ym(2)-ym(1);                % Grid spacing in y
dt=CFL*dx;                     % Timestep based on cfl number (Umax = 1)

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

%a vector for tracking the particles
M = zeros(1250,5);

%% Advancing in time
u_error = 1;
v_error = 1;
limit = 1e-12;

% while u_error > limit && v_error > limit
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
    
    %update u and v star based on above results
    uStar = u;
    vStar = v;
    uStar(2:nx+1,2:ny-1) = uStar(2:nx+1,2:ny-1) + duStar(2:nx+1,2:ny-1);
%     uStar(:,1) = uStar(:,2);
%     uStar(:,end) = uStar(:,end-1);
    vStar(2:nx,2:ny) = vStar(2:nx,2:ny) + dvStar(2:nx,2:ny);
    
    % Make sure, mass is conserved (additively):
    u_diff = u_in - sum(uStar(nx+1,:));
    uStar(nx+1,:) = uStar(nx+1,:)+(u_diff/ny);
    sum(uStar(1,:)-uStar(nx+1,:));
    
    %Solving eqn 17 for del*(del p)
    [ap,ae,aw,an,as,rhs] = eq17(uStar,vStar,dx,dy,dt,nx,ny);
    
    %Set up matrix and use CG to solve for pressure effects
    B=[an ae ap aw as];
    d=[-nx -1 0 1 nx];
    A=(spdiags(B, d, nx*ny, nx*ny))';
    p_cg = solveCG(A,rhs,reshape(p,[nx*ny,1]));
    p_0 = reshape(p_cg,[nx,ny]);
    u_old = u;
    v_old = v;
    [u, v] = eq18(p_0, uStar, vStar,dt,dx,dy);
    
    
    
    % Output the divergence of the velocity field to check solenoidality 
    max_div = max(max((u(2:end,:)-u(1:end-1,:))/dx + ...
        (v(:,2:end)-v(:,1:end-1))/dy));
 
    %update the velocity and position of each particle
     M = particles(M, u, v, i, dx, dy, dt, St);    
    
    %update remaining terms to new values
    Hu_0 = Hu_1;
    Hv_0 = Hv_1;
    t = t + dt
    p = p_0;
    i = i + 1;
    u_error = norm(u - u_old);
    v_error = norm(v - v_old);
    
    %plotting for watching particle flow, slows code dramatically
%     figure(1);
%     clf
%     hold on
%     hs=streamslice(x(2:end),y(2:end),u(1:nx,1:ny)',v(1:nx,1:ny)');
%     set(hs,'color','k','linewidth',1)
%     scatter(M(:,2),M(:,3),20,'filled','r');
%     xlim([0,2]);
%     ylim([0,1]);
    
end

figure(1);
hold on
hs=streamslice(x(2:end),y(2:end),u(1:nx,1:ny)',v(1:nx,1:ny)');
set(hs,'color','k','linewidth',1)
plot(M(:,2),M(:,3),'.','LineWidth',5);
xlim([0,2]);
ylim([0,1]);

figure(2)
plot(ym, u(ceil(nx/2),:))
title('$\mathbf{u}(x=1,y)$ for t=2','fontsize',16,'interpreter','latex') 
ylabel('u','fontsize',12,'interpreter','latex') 
xlabel('y','fontsize',12,'interpreter','latex') 
print('u_1d_t3','-dpng')
