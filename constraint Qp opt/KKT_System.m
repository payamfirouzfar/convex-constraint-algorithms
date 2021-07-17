%% Solving an equality constrained QP problem 

clc, clear, close all




%% QP structure           0.5 * x' * P * x + c' * x      s.t. A * x = b


P = [6 2; 2 5];

c = randn(2,1);

A = [1 1];

b = 3;


[x1,~,~,~,dual_variable1] = quadprog(P,c,[],[],A,b)

[x2, dual_variable2] = KKT_Solve(P, c , A , b )



f = @(x,y) [x y] * P * [x y]' + c' * [x y]';


syms x  y


AffineSet = A * [x y]' - b;



ezcontour(f)

hold on


ezplot(AffineSet)

hold on



plot(x2(1), x2(2), '*')




