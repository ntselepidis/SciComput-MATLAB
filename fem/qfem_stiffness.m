function [W, M, f] = qfem_stiffness(x1, y1, x2, y2, x3, y3, x4, y4)

    y41 = y4 - y1;
    y12 = y1 - y2;
    x14 = x1 - x4;
    x21 = x2 - x1;

    J = 0.25*(x21*y41 - x14*y12);

    W = zeros(4,4);
    W(1,1) = (2*x14^2 + 3*x14*x21 + 2*x21^2 + 2*y12^2 + 3*y12*y41 + 2*y41^2)/(24*J);
    W(2,2) = (2*x14^2 - 3*x14*x21 + 2*x21^2 + 2*y12^2 - 3*y12*y41 + 2*y41^2)/(24*J);
    W(3,3) = W(1,1);
    W(4,4) = W(2,2);
    W(1,2) = -(2*x14^2 - x21^2 - y12^2 + 2*y41^2)/(24*J);
    W(1,3) = -(x14^2 + 3*x14*x21 + x21^2 + y12^2 + 3*y12*y41 + y41^2)/(24*J);
    W(1,4) = (x14^2 - 2*x21^2 - 2*y12^2 + y41^2)/(24*J);
    W(2,3) = W(1,4);
    W(2,4) = -(x14^2 - 3*x14*x21 + x21^2 + y12^2 - 3*y12*y41 + y41^2)/(24*J);
    W(3,4) = W(1,2);
    W(2,1) = W(1,2); W(3,1) = W(1,3); W(4,1) = W(1,4);
    W(3,2) = W(2,3); W(4,2) = W(2,4);
    W(4,3) = W(3,4);

    M = (J/9) * [4 2 1 2; ...
                 2 4 2 1; ...
                 1 2 4 2; ...
                 2 1 2 4];

    f = J * ones(4,1);

end
