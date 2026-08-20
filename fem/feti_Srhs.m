function g = feti_Srhs(bc, FETI)
    Bc = FETI.Bc;
    ndoms = length(Bc);
    nl = size(Bc{1}, 1);
    g = zeros(nl, 1);
    for i = 1 : ndoms
        g = g - Bc{i} * feti_local_solve(i, bc{i}, FETI);
    end
end
