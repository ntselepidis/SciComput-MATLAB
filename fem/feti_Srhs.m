function g = feti_Srhs(bc, FETI)
    [Ac, Bc] = deal( FETI.Ac, FETI.Bc );
    [L, U, P, Q] = deal( FETI.L, FETI.U, FETI.P, FETI.Q );
    ndoms = length(L);
    nl = size(Bc{1}, 1);
    g = zeros(nl, 1);
    for i = 1 : ndoms
        g = g - Bc{i} * sp_solve(L{i}, U{i}, P{i}, Q{i}, bc{i});
        %g = g - Bc{i} * (Ac{i} \ bc{i});
    end
end
