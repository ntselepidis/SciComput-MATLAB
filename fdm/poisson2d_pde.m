clear; clc; close all;
% ------------------------------------- %
%  Delta u(x,y) = f, (x,y) in [a,b]^2   %
% ------------------------------------- %
%  Dirichlet BCs : u(x,y) = 0           %
% ------------------------------------- %

m = 32; % num of intervals per dimension 

a = 0;
b = 1;
h = (b-a)/m;

f = ones((m-1)^2,1); % source term

A1 = (1 / h^2) * five_point(m); % coefficient matrix 5-point molecule

A2 = (1 / h^2) * nine_point(m); % coefficient matrix 9-point molecule

B2 = 6 * speye((m-1)^2) + (1/2) * five_point(m);

u1 = A1\f;

u2 = A2\(B2*f);

figure, mesh( reshape(u1, m-1, m-1) ), title('5-point molecule');
figure, mesh( reshape(u2, m-1, m-1) ), title('9-point molecule');

% figure, spy(A1), title('5-point molecule');
% figure, spy(A2), title('9-point molecule');



