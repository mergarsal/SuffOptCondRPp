function [struct_output] = checkOptimalityCondition(struct_input)
% CHECKOPTIMALITYB4 
% INPUT: 
% - `threshold`:          threshold epsilon in paper
% - `min_sv_relaxation`:  min. singular value of the chosen relaxation
% - `xcc':                eigenvector associated with the smallest eigenvalue of Q
% - `Bs`:                 Wide matrix with all the constraint matrices for this relaxation: 12 x (12 · 5)
% - `Q`:                  data matrix Q (12 x 12)
% - `x`:                  solution to the relaxtive pose problem as [e, t]
% - `f`:                  objective value attain by the solution x
% - 'P':                  scale matrix


% OUTPUT
% - `approx_mu_min`:      approx. of the minimum eigenvalue of M 
% - `is_opt`:             flag: optimality could be certified by bounds
% - 'norm_frob':          norm of the frob matrix A(x_diff)
% - 'norm_Q':             norm (2) of the vector Qprime * x


    
    % TODO
    % check if struct has the desired fields

    Qprime = struct_input.Q; 
    Qprime(10:12, 10:12) = - struct_input.f * eye(3); 

    thresh_prime = struct_input.threshold - struct_input.f;

    delta_x = abs(struct_input.x) - abs(struct_input.xcc);
    X = kron(eye(5), delta_x);
    norm_A_delta = norm(struct_input.Bs * X * struct_input.P, 'fro');


    approx_mu_min = (norm(Qprime * struct_input.x, 2) + thresh_prime * norm_A_delta) / struct_input.min_sv_relaxation; 



    % save to output
    struct_output.approx_mu_min = approx_mu_min;
    struct_output.is_opt = (approx_mu_min < thresh_prime); 
    struct_output.norm_frob = norm_A_delta;
    struct_output.norm_Q = norm(Qprime * struct_input.x, 2); 

end  
