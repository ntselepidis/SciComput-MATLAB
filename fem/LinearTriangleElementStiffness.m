function k = LinearTriangleElementStiffness(E,NU,t,x1,y1,x2,y2,x4,y4,p)

A = (x1*(y2-y4) + x2*(y4-y1) + x4*(y1-y2))/2;

y24 = y2-y4;
y41 = y4-y1;
y12 = y1-y2;
x42 = x4-x2;
x14 = x1-x4;
x21 = x2-x1;

% six dofs per triangle - two dofs at each node
% B is of size 3x6
B = [y24  0  y41  0  y12  0;
      0  x42  0  x14  0  x21;
     x42 y24 x14 y41 x21 y12 ] / (2*A);

if ( p == 1 ) % plane stress
    
    D = (E/(1-NU*NU))*[1 NU 0; NU 1 0; 0 0 (1-NU)/2];

elseif ( p == 2 ) % plane strain

    D = (E/(1+NU)/(1-2*NU))*[1-NU NU 0; NU 1-NU 0; 0 0 (1-2*NU)/2];

end

k = t*A*B'*D*B;

end