function g = Srhs(b, SC)

    [In, out, Aoi] = deal( SC.In, SC.out, SC.Aoi );
    [L, U, P, Q] = deal( SC.L, SC.U, SC.P, SC.Q );    
    ndoms = length(L);

    g = b(out);
    for i=1:ndoms
        g = g - Aoi{i} * sp_solve( L{i}, U{i}, P{i}, Q{i}, b(In{i}) ); 
    end
    
end