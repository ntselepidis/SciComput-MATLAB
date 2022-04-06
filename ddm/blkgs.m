function xP = blkgs(xP, params, AP, bP)

    blk   = params.blk;
    L     = params.L;
    U     = params.U;
    P     = params.P;
    Q     = params.Q;
    ndoms = length(L);
    for i = 1 : ndoms
        rP = bP - AP*xP;
        %omega = (rP'*rP) / (rP'*AP*rP);
        omega = 1.7;
        xP( blk(i):blk(i+1)-1 ) = xP( blk(i):blk(i+1)-1 ) ...
            + omega * sp_solve(L{i}, U{i}, P{i}, Q{i}, rP( blk(i):blk(i+1)-1 ));
    end

end
