function y = feti_Smultx(x, FETI)
    Bc = FETI.Bc;
    ndoms = length(Bc);
    nl = size(Bc{1}, 1);
    y = zeros(nl, 1);
    for i = 1 : ndoms
        y = y - Bc{i} * feti_local_solve(i, Bc{i}'*x, FETI);
    end
end
