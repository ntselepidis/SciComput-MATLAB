clear; clc; close all;
rng(0); % used for reproducibility
[a, b] = deal(-1, 0);
m = 64; % number of subintervals per dimension (per square domain)
fem_type = 0; % 0 for triangle and 1 for square
if ( fem_type )
    [h, ne, n, coo1, con1, bounds1] = qfem_discretize(a, b, m);
else
    [h, ne, n, coo1, con1, bounds1] = tfem_discretize(a, b, m);
end
[coo, con, bounds, ~, ~] = get_lshape_from_square(a, b, m, coo1, con1);

f = 0; % Helmholtz -> Poisson
if (fem_type)
    [A, b] = qfem_assemble(coo, con, f);
else
    [A, b] = tfem_assemble(coo, con, f);
end
disp('Matrix Assembly - Done');

% Domain Decomposition
ndoms = 10;
map = kmeans(coo, ndoms);
disp('kmeans - Done');
clr = ['b' 'r' 'g' 'm' 'c'];
while(length(clr) < ndoms)
    clr = [clr clr];
end
In = cell(ndoms,1);
Out = cell(ndoms,1);
G = abs(A)+abs(A)';
figure, hold on,
for i = 1:ndoms
    dpnts = find(map == i);
    plot(coo(dpnts,1), coo(dpnts,2), strcat(clr(i), '*'));
    for jj = 1:length(dpnts)
        j = dpnts(jj);
        nbr = setdiff(find(G(:,j)), j);
        if (any(map(nbr) ~= i))
            Out{i} = [Out{i} j];
        else
            In{i} = [In{i} j];
        end
    end
    In{i} = setdiff(In{i}, bounds);
    Out{i} = setdiff(Out{i}, bounds);
end
hold off;

in = horzcat(In{:});
out = horzcat(Out{:});
p = [in out];
disp('Data Partitioning - Done');

% Solve by Schur Complement Method
xx = zeros(n,1);
% S = A(out,out)-A(out,in)*(A(in,in)\A(in,out));
% g = b(out)-A(out,in)*(A(in,in)\b(in));
S = A(out,out);
g = b(out);
for i = 1:ndoms
    S = S-A(out,In{i})*(A(In{i},In{i})\A(In{i},out));
    g = g-A(out,In{i})*(A(In{i},In{i})\b(In{i}));
end
% xout = S\g;
[L,U] = ilu(S);
lfill = 1; 
M = make_spai(S, lfill);
xout = bicgstab(S, g, 1e-10, 500);
xout = bicgstab(S, g, 1e-10, 500, @(y) U\(L\y));
xout = bicgstab(S, g, 1e-10, 500, @(y) M*y);
xout = gmres(S, g, 20, 1e-10, 100);
xout = gmres(S, g, 20, 1e-10, 100, @(y) U\(L\y));
xout = gmres(S, g, 20, 1e-10, 100, @(y) M*y);
% xin = A(in,in)\(b(in)-A(in,out)*xout);
for i = 1:ndoms
    xx(In{i}) = A(In{i},In{i})\(b(In{i})-A(In{i},out)*xout);
end
% xx(in) = xin;
xx(out) = xout;
xx(bounds) = 0;
figure, trimesh(con, coo(:,1), coo(:,2), xx);
