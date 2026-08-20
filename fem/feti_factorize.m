function [L, U, P, Q] = feti_factorize(Ac)
    
    ndoms = length(Ac);
    L = cell(ndoms, 1);
    U = cell(ndoms, 1);
    P = cell(ndoms, 1);
    Q = cell(ndoms, 1);
    for i = 1 : ndoms
        [L{i}, U{i}, P{i}, Q{i}] = lu( Ac{i}, 'vector' );
    end

end
