function SC = sc_prep(A, In, Out)
ndoms = length(In);
out = horzcat(Out{:});
Aoo = A(out,out); 
[Aio, Aoi, L, U, P, Q] = deal( cell(ndoms,1) );
for i = 1 : ndoms
    inn = In{i};
    Aio{i} = A(inn,out);
    Aoi{i} = A(out,inn);
    [L{i}, U{i}, P{i}, Q{i}] = lu( A(inn,inn), 'vector' );
end
[SC.L, SC.U, SC.P, SC.Q] = deal(L, U, P, Q);
[SC.Aoo, SC.Aio, SC.Aoi] = deal(Aoo, Aio, Aoi);
[SC.In, SC.Out, SC.out] = deal(In, Out, out);
end
