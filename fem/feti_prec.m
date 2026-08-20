function y = feti_prec(x, FETI)
    [Ac, Bc, Bdc, Sc, Bdco] = deal( FETI.Ac, FETI.Bc, FETI.Bdc, FETI.Sc, FETI.Bdco );
    ndoms = length(Ac);
    nl = size(Bc{1}, 1);
    y = zeros(nl, 1);
    for i = 1 : ndoms
        %y = y - Bdc{i} * (Ac{i} * (Bdc{i}'*x));
        y = y - Bdco{i} * (Sc{i} * (Bdco{i}'*x));
    end
end
