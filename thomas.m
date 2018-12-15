function y = thomas(aw,ap,ae,b,args,nx)
% Solve a tridiagonal linear system using the Thomas Algorithm.

v = zeros(args,1);
y = v;
w = ap(1);
y(1) = b(1)/w;
for i = 2:args
    v(i-1) = ae(i-1)/w;
    w = ap(i) - aw(i)*v(i-1);
    y(i) = (b(i) - aw(i)*y(i-1))/w;         
end
for j = args-1:-1:1
    y(j) = y(j) - v(j)*y(j+1);
end
end