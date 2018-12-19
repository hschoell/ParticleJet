% Just to test the SOR solver
aw = [0 3 3 3 0 3 3 3 0 3 3 3];
ae = [3 3 3 0 3 3 3 0 3 3 3 0];
an = [3 3 3 3 3 3 3 3 0 0 0 0];
as = [0 0 0 0 3 3 3 3 3 3 3 3];
ap = [5 5 5 5 5 5 5 5 5 5 5 5];

B=[an' ae' ap' aw' as'];
d=[-4 -1 0 1 4];
A=(spdiags(B, d, 12, 12))';
A1 = full(A)
x = [6 8 4 5 90 5 3 76 9 65 7 9]'
rhs = A*x;
x_0 = A\rhs

[x_1,res] = solveSOR(aw,ae,an,as,ap,rhs,zeros(4,3),4,3)