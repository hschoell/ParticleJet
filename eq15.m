function duStar = eq15(duStarStar, d2, dt, Re)
% This function calculates duStar based uporn equation 15

% Determine Number of unknowns
args = length(duStarStar(1,:))*length(duStarStar(:,1));
nx = length(duStarStar(:,1));

% Set up the operator on the left-hand-side
as = -dt/(2*d2^2*Re)*ones(args,1);
ap = (1+dt/(d2^2*Re))*ones(args,1);
an = -dt/(2*d2^2*Re)*ones(args,1);

% Set up the right-hand-side. Note, that now inversed lexiographical
% ordering is used. (Vector should look like [u(1,1) u(2,1) u(3,1) ...]
rhs = reshape(duStarStar',[1,args])';

% Solve the system using the Thomas-Algorithm
duStar_vec = thomas(as,ap,an,rhs,args,nx);
duStar = reshape(duStar_vec,[length(duStarStar(1,:)),length(duStarStar(:,1))])';
end
