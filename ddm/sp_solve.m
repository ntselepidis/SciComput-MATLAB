function x = sp_solve(L, U, p, q, b)
x = U \ (L \ b(p, :));
x(q, :) = x;
end
