%% Problem Structure

%           min f(x) = - \sum_{i = 1} ^ {n} log(x_i)
%       
%           s.t. A * x = b

clc, clear, close all

A = [1 1];

b = 6;

eps = 1e-5;

x = [0.1 5.9]';


options = struct('Maxiter', 1000, 'tolerance', eps, 'Initial_Condition', x);


[xopt, fval, Iter, X ] = EQ_NM(A, b, options );


sdpvar h(2,1)

obj = -log(h(1)) - log(h(2));

const = A * h == b;

a = optimize(const,obj);

h = value(h);





f = @(x,y) -log(x) - log(y); 



syms x1 x2

AffineSet = A * [x1 x2]' - b;



ezcontour(f)
hold on


ezplot(AffineSet)
hold on


plot(x(1),x(2),'o')
hold on


plot(xopt(1) , xopt(2),'*')
hold on


plot(X(1,:), X(2,:),'r*-')


legend('contour of f(x)', 'feasible set', 'initial condition','optimal solution','iterations')













