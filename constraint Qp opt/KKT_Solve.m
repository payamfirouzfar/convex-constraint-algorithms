%% Solving KKT system by block elimination method

%       [P      A';A        0] * [x     v]' = [-c       b]
%
%       P * x + A' * v = -c ===> x = - inv(P) * (A' * v + c)

function [Primal_Optimal, Dual_Optimal] = KKT_Solve(P, c , A , b )


inv_P = inv(P);

S = A * inv_P * A';


Dual_Optimal = - S \ (A * inv_P * c + b);


Primal_Optimal = - P \(A' * Dual_Optimal + c);



end

