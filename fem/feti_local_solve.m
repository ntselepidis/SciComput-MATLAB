function y = feti_local_solve(i, b, FETI)
    [L, U, P, Q, R] = deal( FETI.L, FETI.U, FETI.P, FETI.Q, FETI.R );
    if ~FETI.floating(i)
        y = sp_solve(L{i}, U{i}, P{i}, Q{i}, b);
    else
        rhs = [b; zeros(size(R{i},2),1)];
        z = sp_solve(L{i}, U{i}, P{i}, Q{i}, rhs);
        y = z(1:length(b));
    end
end

