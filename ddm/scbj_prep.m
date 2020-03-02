function SBJ = scbj_prep(A,SC)
[In,Out,L,U,P,Q] = deal(SC.In,SC.Out,SC.L,SC.U,SC.P,SC.Q);
ndoms = length(In);
[S, Ls, Us, Ps, Qs, T] = deal( cell(ndoms,1) );
blk2 = zeros(ndoms+1,1);
blk2(1) = 1;
for i=1:ndoms
%     S{i} = A(Out{i},Out{i}) - A(Out{i},In{i}) * ...
%         sp_solve( L{i},U{i},P{i},Q{i}, A(In{i},Out{i}) );
    T{i} = A(Out{i},In{i}) * sp_solve( L{i},U{i},P{i},Q{i}, A(In{i},Out{i}) );
    S{i} = A(Out{i},Out{i}) - T{i};
    [Ls{i}, Us{i}, Ps{i}, Qs{i}] = lu( S{i}, 'vector' );
    blk2(i+1) = length(Out{i});
end
blk2 = cumsum(blk2);
[SBJ.L, SBJ.U, SBJ.P, SBJ.Q, SBJ.blk] = deal( Ls, Us, Ps, Qs, blk2 );
SBJ.T = T;
end