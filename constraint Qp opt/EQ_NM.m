%% Problem Structure

%           min f(x) = - \sum_{i = 1} ^ {n} log(x_i)
%       
%           s.t. A * x = b
%
%           drev{f(x)} = [-1/x_1 ,...., -1/x_n]'
%
%           Hessian = diag(1/(x_1^2),...,1/(x_n^2))

function [xopt, fval, Iter, X ] = EQ_NM(A, b, options )

%% Algorithm Parameters

MaxIter = options.Maxiter;


x = options.Initial_Condition;


eps = options.tolerance;


t = 0.2;


%% Main Loop


for k = 1: MaxIter
    
    
    Gradient = [-1/x(1) -1/x(2)]';
    
    Hessian = blkdiag(1/x(1)^2 , 1/x(2)^2);
    
    
    [v, w] = KKT_Solve(Hessian, Gradient , A , zeros(size(A,1),1));
    
    
    X(:,k) = x;
    
    x = x + t * v;
    
    r = A * x - b;
    
    Drev_L = Gradient + A' * w;
    
    if norm(r) <= eps && norm(Drev_L)<= eps
        
        
        xopt = x;
        
        fval = - (log(x(1)) + log(x(2)));
        
        Iter = k;
        
        break
        
    elseif k == MaxIter
        
        error('Please increase number of iterations')
        
    end
    
    
    
end












end

