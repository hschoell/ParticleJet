% Function to calculate del(del p) using 
function [ap,ae,aw,an,as,rhs] = eq17(uStar,vStar,u,v,dx,dy,dt,nx,ny);

% Number of unknowns
args = length(u(1,:))*length(v(:,1));

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

for m=1:128
    for n = 1:64
        d = m+128*(n-1);
        if m == 1 &&  n ~= 1
            aw(d) = 0;
        elseif m == 128 && n ~=64
            ae(d) = 0;
        end
        if n == 1 && m ~= 128
            as(d) = 0;
        elseif n == 64 && m ~= 128
            an(d) = 0;
        end
    end
end

for s=1:128
    for t = 1:64
        e = s+128*(t-1);
        aw1(s,t) = aw(e);
        as1(s,t) = as(e);
        ae1(s,t) = ae(e);
        an1(s,t) = an(e);
    end
end

end







