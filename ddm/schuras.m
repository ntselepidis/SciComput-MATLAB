function x = schuras(b, params)

    overlap = params.overlap; % masked in preprocessing
    L       = params.L;
    U       = params.U;
    P       = params.P;
    Q       = params.Q;
    ndoms   = length(L);
    
    x = zeros( length(b), 1 );
    for i = 1 : ndoms
        bb = b( overlap{i} );
        xx = sp_solve( L{i}, U{i}, P{i}, Q{i}, bb );
        x( overlap{i} ) = x( overlap{i} ) + xx;
    end
    
end
