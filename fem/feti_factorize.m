function [L, U, P, Q] = feti_factorize(Ac, FETI)

ndoms = length(Ac);

L = cell(ndoms,1);
U = cell(ndoms,1);
P = cell(ndoms,1);
Q = cell(ndoms,1);
R = FETI.R;

for i = 1:ndoms

    if ~FETI.floating(i)

        % Non-floating domain
        [L{i},U{i},P{i},Q{i}] = lu(Ac{i},'vector');

    else

        % Floating domain:
        % [ A R ; R' 0 ]
        r = size(R{i},2);

        K = [Ac{i}, R{i};
             R{i}', sparse(r,r)];

        [L{i},U{i},P{i},Q{i}] = lu(K,'vector');

    end
end
end
