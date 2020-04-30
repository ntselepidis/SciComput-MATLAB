function x = mgmu3(A, b, x, P, R, M, n1, n2, lev, levmax, omega, mu)
    if ( lev == levmax )
        x{lev} = A{lev} \ b{lev};
        return;
    end
    for i = 1:n1
        x{lev} = x{lev}+omega*(M{lev}\(b{lev}-A{lev}*x{lev}));
    end
    b{lev+1} = R{lev}*(b{lev}-A{lev}*x{lev});
    x{lev+1} = zeros(length(A{lev+1}), 1);
    for i = 1:mu
        x = mgmu3(A, b, x, P, R, M, n1, n2, lev+1, levmax, omega, mu);
    end
    x{lev} = x{lev}+P{lev}*x{lev+1};
    for i = 1:n2
        x{lev} = x{lev}+omega*(M{lev}\(b{lev}-A{lev}*x{lev}));
    end
end
