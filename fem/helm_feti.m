clear; clc; close all; %rng(1);
[a, b] = deal(0, 1);

f = 0*5e8;    % if f = 0 then Helmholtz -> Poisson

m = 64;       % number of subintervals per dimension
ndoms = 16;    % number of subdomains
fem_type = 0; % 0 -> triangular elements, 1 -> square elements
lshape = 1;   % 0 -> square domain, 1-> L-shape domain
lumped = 0;

% discretize domain
if ( fem_type )
    [h, ne, n, coo, con, bounds] = qfem_discretize(a, b, m);
else
    [h, ne, n, coo, con, bounds] = tfem_discretize(a, b, m);
end

if ( lshape )
    [coo, con, bounds, ~, ~] = get_lshape_from_square(a, b, m, coo, con);
    ne = size(con,1);
end

% compute element centers and partition domain geometrically
ce = zeros(ne,2);
for i = 1:ne
    xc = sum( coo( con(i,:), 1 ) ) / size(con,2);
    yc = sum( coo( con(i,:), 2 ) ) / size(con,2);
    ce(i,:) = [xc yc];
end
part = kmeans(ce, ndoms); % domain partitioning

% plot partitioning
clr = ['b' 'r' 'g' 'c' 'm' 'y']';
clr = repmat(clr, ceil(ndoms/length(clr)), 1);
figure, hold on;
for i = 1:ndoms
    ind = (part==i);
    plot(ce(ind,1), ce(ind,2), strcat(clr(i), '*'));
end
hold off;

elems = (1:ne)';        % element indices
ielem = cell(ndoms,1);  % element indices of each domain
nelem = zeros(ndoms,1); % element number per domain
dcon1 = cell(ndoms,1);  % domain connectivity matrix    (global numbering)
dcon2 = cell(ndoms,1);  % domain connectivity matrix    ( local numbering)
All = cell(ndoms,1);    % vertex indices of each domain (global numbering) 
Ac = cell(ndoms,1);     % subdomain stiffness matrices 
bc = cell(ndoms,1);     % subdomain rhs vectors
xc = cell(ndoms,1);     % subdomain solution vectors
blk = zeros(ndoms+1,1); % subdomain index bounds
blk(1) = 1;
has_dirichlet = false(ndoms,1);

% assemble subdomain stiffness matrices
for i = 1:ndoms
    ielem{i} = elems(part==i);
    nelem(i) = length(ielem{i});
    dcon1{i} = con(ielem{i},:);
    
    All{i} = unique( dcon1{i} );
    mapp = zeros( length(coo), 1 );
    mapp(All{i}) = 1:length(All{i});
    dcon2{i} = mapp( dcon1{i} );
    
    nv = length(All{i});
    dcoo = zeros(nv,2);
    for j = 1:nv
        dcoo(j,:) = coo(All{i}(j),:);
    end
    
    if ( fem_type )
        [Ac{i}, bc{i}] = qfem_assemble(dcoo, dcon2{i}, f);
    else
        [Ac{i}, bc{i}] = tfem_assemble(dcoo, dcon2{i}, f);
    end
    
    % remove boundaries ( dirichlet (0) )
    dbnds = mapp(bounds);
    dbnds = dbnds(dbnds>0);
    has_dirichlet(i) = ~isempty(dbnds); % for determining floating subdomains
    Ac{i}(dbnds,:) = [];
    Ac{i}(:,dbnds) = [];
    bc{i}(dbnds) = [];
    All{i} = setdiff(All{i}, bounds);
    blk(i+1) = blk(i)+length(All{i});
end

p = vertcat(All{:});
n = size(coo,1);
deg = zeros(n,1);
for i = 1:ndoms
    deg(All{i}) = deg(All{i})+1;
end
s = find(deg>1);
deg = deg(deg>1);

% serial B and Bd computation
B = spalloc(sum(deg-1), length(p), sum( 2*(deg-1) ));
Bd = spalloc(sum(deg-1), length(p), sum( deg.*(deg-1) ));
idx = 1;
for i = 1:length(s)
    ind = find(p==s(i));
    count = 1;
    for j = 1:deg(i)-1
        B(idx,ind(j)) = 1;
        B(idx,ind(j+1)) = -1;
        for k = 1 : j
            Bd(idx,ind(k)) = 1-count/deg(i);
        end
        for k = j+1 : length(ind)
            Bd(idx,ind(k)) = -count/deg(i);
        end
        count = count+1;
        idx = idx+1;
    end
end

Bc = cell(ndoms,1);
Bdc = cell(ndoms,1);
for i = 1:ndoms
    idx = blk(i):blk(i+1)-1;
    Bc{i} = B(:,idx);
    Bdc{i} = Bd(:,idx);
end

%nl = size(B,1);
%S = zeros(nl,nl);
%g = zeros(nl,1);
%for i = 1:ndoms
%    S = S + Bc{i} * ( Ac{i} \ Bc{i}' );
%    g = g + Bc{i} * ( Ac{i} \ bc{i}  );
%end

% AA = blkdiag(Ac{:});
% lambda = pcg(-S, -g, 1e-10, 100);
% lambda = pcg(-S, -g, 1e-10, 100, @(y) -B*(AA*(B'*y)));
% lambda = S\g;
% lambda = 0*lambda;

% Reorder based on In and Out
Lin = cell(ndoms,1);
Lout = cell(ndoms,1);
reord = cell(ndoms,1);
Sc = cell(ndoms,1);
Bdco = cell(ndoms,1);
for i = 1 : ndoms
    % Find local in and out points
    [~, Lout{i}] = find(Bc{i});
    Lout{i} = unique(Lout{i});
    Lin{i} = setdiff((1:length(Ac{i}))', Lout{i});

    % Compute local Schur complements and slice Bdc
    if ~lumped
        Sc{i} = Ac{i}(Lout{i}, Lout{i}) - Ac{i}(Lout{i}, Lin{i}) * ( Ac{i}(Lin{i}, Lin{i}) \ Ac{i}(Lin{i}, Lout{i}) );
    else
        Sc{i} = Ac{i}(Lout{i}, Lout{i});
    end
    Bdco{i} = Bdc{i}(:, Lout{i});

    % Reorder data in place
    reord{i} = [Lin{i}; Lout{i}];
    All{i} = All{i}(reord{i});
    Ac{i} = Ac{i}(reord{i}, reord{i});
    Bc{i} = Bc{i}(:, reord{i});
    Bdc{i} = Bdc{i}(:, reord{i});
    bc{i} = bc{i}(reord{i});
end

% Floating subdomains: scalar Poisson rigid-body modes
R = cell(ndoms,1);
for i = 1:ndoms
    if ~has_dirichlet(i)
        R{i} = ones(length(All{i}),1);
        R{i} = R{i}/norm(R{i});
    else
        R{i} = zeros(length(All{i}),0);
    end
end

% FETI coarse/nullspace matrix G = B R
G = [];
for i = 1:ndoms
    if ~isempty(R{i})
        G = [G, Bc{i}*R{i}];
    end
end

FETI.floating = ~has_dirichlet;
FETI.R = R;
FETI.G = G;

% Factorize domains
[L, U, P, Q] = feti_factorize(Ac, FETI);
[FETI.L, FETI.U, FETI.P, FETI.Q] = deal(L, U, P, Q);

% Setup FETI struct
[FETI.Ac, FETI.Bc, FETI.Bdc, FETI.Sc, FETI.Bdco] = deal(Ac, Bc, Bdc, Sc, Bdco);

Z = @(x) feti_project(x, FETI);
projected_matvec = @(x) Z(feti_Smultx(Z(x), FETI));
projected_precon = @(x) Z(feti_prec(x, FETI));
projected_rhs = Z(feti_Srhs(bc, FETI));

% Run PCG
%lambda = pcg(@(x) feti_Smultx(x, FETI), feti_Srhs(bc, FETI), 1e-10, 100);
%lambda = pcg(@(x) feti_Smultx(x, FETI), feti_Srhs(bc, FETI), 1e-10, 100, @(x) feti_prec(x, FETI));
lambda = pcg(@(x) projected_matvec(x), projected_rhs, 1e-10, 100, @(x) projected_precon(x));

% Backsolve
for i = 1:ndoms
    %xc{i} = Ac{i} \ ( bc{i} - Bc{i}'*lambda ); 
    %xc{i} = sp_solve(L{i}, U{i}, P{i}, Q{i}, bc{i} - Bc{i}'*lambda); 
    xc{i} = feti_local_solve(i, bc{i} - Bc{i}'*lambda, FETI); 
end

if ~isempty(G)
    % Fix floating subdomain jumps
    % 1. Compute interface mismatch
    jump = zeros(size(lambda));
    for i = 1:ndoms
        jump = jump + Bc{i}*xc{i};
    end

    % 2. Find best alpha
    alpha = -(G'*G)\(G'*jump);

    % 3. Add constants to shift floating subdomains
    col = 0;
    for i = 1:ndoms
        if FETI.floating(i)
            col = col + 1;
            xc{i} = xc{i} + R{i}*alpha(col);
        end
    end
end

x = zeros(n,1);
for i = 1:ndoms
    x(All{i}) = xc{i};
end

figure, trimesh(con, coo(:,1), coo(:,2), x), % plot solution
title(sprintf('Solution of Helmholtz PDE for f = %e',f));

xx=vertcat(xc{:});
BB=horzcat(Bc{:});
norm(BB*xx)
