function [L, U, P, Q] = factorize(A, All)
    
    ndoms = length(All);
    L = cell(ndoms, 1);
    U = cell(ndoms, 1);
    P = cell(ndoms, 1);
    Q = cell(ndoms, 1);
    for i = 1 : ndoms
        [L{i}, U{i}, P{i}, Q{i}] = lu( A(All{i},All{i}), 'vector' );
    end

end
