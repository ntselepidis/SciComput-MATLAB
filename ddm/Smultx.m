function y = Smultx(x, SC)

    [Aoo, Aio, Aoi] = deal(SC.Aoo, SC.Aio, SC.Aoi);
    [L, U, P, Q] = deal( SC.L, SC.U, SC.P, SC.Q );
    ndoms = length(L);

    y = Aoo*x;
    for i=1:ndoms
        y = y - Aoi{i} * sp_solve( L{i}, U{i}, P{i}, Q{i}, Aio{i} * x ); 
    end
    
end