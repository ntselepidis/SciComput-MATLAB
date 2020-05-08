function x = blkjac(b, params)

    blk   = params.blk;
    L     = params.L;
    U     = params.U;
    P     = params.P;
    Q     = params.Q;
    ndoms = length(L);
    bc    = cell(ndoms, 1);
    xc    = cell(ndoms, 1);
    for i = 1 : ndoms
        bc{i} = b( blk(i):blk(i+1)-1 );
        xc{i} = sp_solve(L{i}, U{i}, P{i}, Q{i}, bc{i});
    end
    x = vertcat(xc{:});

end
