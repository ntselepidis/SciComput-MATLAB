function [W, M, f] = tfem_stiffness(x1, y1, x2, y2, x4, y4)
    y41 = y4 - y1;
    y12 = y1 - y2;
    x14 = x1 - x4;
    x21 = x2 - x1;

    J = 0.25*((x2 - x1)*(y4 - y1) - (x4 - x1)*(y2 - y1));

    W = zeros(3,3);
    W(1,1) = (x14^2 + 2*x14*x21 + x21^2 + y12^2 + 2*y12*y41 + y41^2)/(8*J);
    W(2,2) = (x14^2 + y41^2)/(8*J);
    W(3,3) = (x21^2 + y12^2)/(8*J);
    W(1,2) = -(x14^2 + x21*x14 + y41^2 + y12*y41)/(8*J);
    W(1,3) = -(x21^2 + x14*x21 + y12^2 + y41*y12)/(8*J);
    W(2,3) = (x14*x21 + y12*y41)/(8*J);
    W(2,1) = W(1,2); W(3,1) = W(1,3);
    W(3,2) = W(2,3);

    M = zeros(3,3);
    M(1,1) = J/3;
    M(2,2) = J/3;
    M(3,3) = J/3;
    M(1,2) = J/6;
    M(1,3) = J/6;
    M(2,3) = J/6;
    M(2,1) = M(1,2); M(3,1) = M(1,3);
    M(3,2) = M(2,3);

    f = ((2*J)/3)*ones(3,1);
end