function x = mti_prec(b, MTI)

[All, Bc, Bdc] = deal( MTI.All, MTI.Bc, MTI.Bdc );
[L, U, P, Q] = deal( MTI.L, MTI.U, MTI.P, MTI.Q );

ndoms = length(All);
nlm = size(Bc{1},1);
lm  = zeros(nlm,1);
xc  = cell(ndoms,1);
x   = zeros(length(b),1);
for i = 1 : ndoms
    xc{i} = sp_solve( L{i},U{i},P{i},Q{i}, b(All{i}) );
    lm = lm + Bdc{i}*xc{i};
end
for i = 1 : ndoms
    xc{i} = xc{i} - Bc{i}'*lm;
    x(All{i}) = xc{i};
end

end
