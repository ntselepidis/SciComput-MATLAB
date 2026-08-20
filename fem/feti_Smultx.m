function y = feti_Smultx(x, FETI)
    [Ac, Bc] = deal( FETI.Ac, FETI.Bc );
    [L, U, P, Q] = deal( FETI.L, FETI.U, FETI.P, FETI.Q );
    ndoms = length(L);
    nl = size(Bc{1}, 1);
    y = zeros(nl, 1);
    for i = 1 : ndoms
        y = y - Bc{i} * sp_solve(L{i}, U{i}, P{i}, Q{i}, Bc{i}'*x);
        %y = y - Bc{i} * (Ac{i} \ (Bc{i}'*x));
    end
end
