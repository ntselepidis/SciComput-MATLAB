clear; clc; close all;
[a, b] = deal(-1, 0);
m = 25; % number of subintervals per dimension (per square domain)
fem_type = 1; % 0 for triangle and 1 for square
if ( fem_type )
    [h, ne, n, coo1, con1, bounds1] = qfem(a, b, m);
else
    [h, ne, n, coo1, con1, bounds1] = tfem(a, b, m);
end
[coo, con, bounds, sep, dpnts] = get_lshape_from_square(a, b, m, coo1, con1);

figure, 
plot(coo(:,1), coo(:,2), '*', ...
coo(sep,1), coo(sep,2), 'go', ...
coo(bounds,1), coo(bounds,2), 'r*');

f = 0; % Helmholtz -> Poisson
[A, b] = qfem_assemble(coo, con, f);

% Domain Decomposition
% use three subdomains (fixed)
In = cell(3,1);
Out = cell(3,1);
sep = setdiff(sep, bounds);
for i = 1:3
    dpnts{i} = setdiff(dpnts{i}, bounds);
    In{i} = setdiff(dpnts{i}, sep);
    Out{i} = setdiff(dpnts{i}, In{i});
end
in = horzcat(In{:});
out = unique(horzcat(Out{:}));
p = [in out];

% Solve by Schur Complement Method
S = A(out,out)-A(out,in)*(A(in,in)\A(in,out));
g = b(out)-A(out,in)*(A(in,in)\b(in));
xout = S\g;
xin = A(in,in)\(b(in)-A(in,out)*xout);

xx = zeros(n,1);
xx(in) = xin;
xx(out) = xout;
xx(bounds) = 0;
figure, trimesh(con, coo(:,1), coo(:,2), xx);
