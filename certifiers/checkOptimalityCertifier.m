function str_out = checkOptimalityCertifier(str_in)


e_min_svd = str_in.e_min_svd;
Q = str_in.Q; 
E = str_in.E; 
t = str_in.t; 
f_hat = str_in.f_hat;  
e_opt = E(:); 
As = str_in.As;   
    
    isOpt=0;
    
        % 2.1 COmpute Lagrange multipliers
        % implement reduce version
        
        Qf = Q - f_hat * blkdiag(zeros(9), eye(3)); 
        
        [L_hat, A_x, b, res]=computeLagrangeMultipliersReducedRelaxation(Qf, As, e_opt, t);
        
        
        % 2.3 Compute min eigenvalu of M
        M = Qf;
        for i=1:size(L_hat, 1)
            M = M - L_hat(i)*As(:, :, i);
        end %end-for each constraint
        d=eig(M);
        mu_min=double(real(min(d)));
        
        % 2.5 Check conditions for optimality
        if double(mu_min > e_min_svd)
            isOpt=1;
        end
        
        str_out = struct(); 
        str_out.is_opt = isOpt; 
        str_out.mu_min = mu_min; 
        str_out.eig_M = d; 
        str_out.mult_hat = L_hat; 
        str_out.res_ls = res;
        str_out.Hessian = M;

   
end

