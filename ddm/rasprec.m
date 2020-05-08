function x = rasprec(b, params)
    
    L     = params.L;
    U     = params.U;
    P     = params.P;
    Q     = params.Q;
    All   = params.All;
    All2  = params.All2;
    ndoms = length(L);
    
    x = zeros(length(b),1);
    for i = 1 : ndoms
        bb = b( All2{i} );
        xx = sp_solve( L{i}, U{i}, P{i}, Q{i}, bb );
        x( All{i} ) = xx( 1:length(All{i}) );
    end
    
end
