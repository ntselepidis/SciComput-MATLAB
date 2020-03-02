function SBP = scbp_prep(A,SC,prob)
[In,Out,L,U,P,Q] = deal(SC.In,SC.Out,SC.L,SC.U,SC.P,SC.Q);
ndoms = length(In);
[Sp, Lp, Up, Pp, Qp] = deal( cell(ndoms,1) );
blk2 = zeros(ndoms+1,1);
blk2(1) = 1;
for i=1:ndoms
    [AOO,AOI,AIO] = deal( A(Out{i},Out{i}), A(Out{i},In{i}), A(In{i},Out{i}) );
    [LL,UU,PP,QQ] = deal( L{i}, U{i}, P{i}, Q{i} );
    Sx = @(x) AOO*x-AOI*sp_solve(LL,UU,PP,QQ,AIO*x);
    Sp{i} = probe2(Sx, prob, length(Out{i}));
    [Lp{i}, Up{i}, Pp{i}, Qp{i}] = lu( Sp{i}, 'vector' );
    blk2(i+1) = length(Out{i});
end
blk2 = cumsum(blk2);
[SBP.L, SBP.U, SBP.P, SBP.Q, SBP.blk] = deal( Lp, Up, Pp, Qp, blk2 );
end