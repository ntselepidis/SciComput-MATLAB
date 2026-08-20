function y = feti_prec(x, FETI)
    [Sc, Bdco] = deal( FETI.Sc, FETI.Bdco );
    ndoms = length(Bdco);
    nl = size(Bdco{1}, 1);
    y = zeros(nl, 1);
    for i = 1 : ndoms
        y = y - Bdco{i} * (Sc{i} * (Bdco{i}'*x));
    end
end
