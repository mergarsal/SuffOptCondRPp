function struct_output = computeMinSingularValueB4All(xcc, Ps, discard_column)
        
        % display('Computing min SV with B4'); 
        
        % check other relaxations 
        x1=xcc(1);x2=xcc(2);x3=xcc(3);x4=xcc(4);x5=xcc(5);x6=xcc(6);
        x7=xcc(7);x8=xcc(8);x9=xcc(9);x10=xcc(10);x11=xcc(11);x12=xcc(12);


        % create matrix A wint coeff.
        A = zeros(12, 7);

        % to simplify
        L1=1;L2=2;L3=3;L4=4;L5=5;L6=6;L7=7;

        % first equation
        A(1, L2) = x1;   A(1, L3) = x2*0.5;  A(1, L4) = x3*0.5;

        % Second equation
        A(2, L3) = x1*0.5;  A(2, L5)=x2; A(2, L6)=x3*0.5;

        % 3rd equation
        A(3, L4) = x1*0.5;  A(3, L6)=x2*0.5; A(3, L7) = x3;

        % 4th equation
        A(4, L2) = x4; A(4, L3)=x5*0.5; A(4, L4) = x6*0.5;

        % 5th equation
        A(5, L3)=x4*0.5; A(5, L5) = x5;   A(5, L6) = x6*0.5;

        % 6th equation
        A(6, L4)=x4*0.5;  A(6, L6) = x5*0.5;  A(6, L7) = x6;

        % 7th equation
        A(7, L2)=x7;  A(7, L3)=x8*0.5;  A(7, L4)=x9*0.5;

        % 8th equation
        A(8, L3)=x7*0.5;  A(8, L5)=x8; A(8, L6)=x9*0.5;

        % 9th equation
        A(9, L4)=x7*0.5;  A(9, L6)=x8*0.5;  A(9, L7)=x9;

        % 10th equation
        A(10, L3)=x11*0.5;  A(10, L4)=x12*0.5;  A(10, L5)=-x10; A(10, L1)=x10; A(10, L7)=-x10;

        % 11th equation
        A(11, L3)=x10*0.5;  A(11, L6)=x12*0.5;  A(11, L2)=-x11; A(11, L1)=x11; A(11, L7)=-x11;

        % 12th equation
        A(12, L4)=x10*0.5;  A(12, L6)=x11*0.5;  A(12, L2)=-x12; A(12, L1)=x12; A(12, L5)=-x12;

        % find the best one 
        best_min_sv = -1000000; 
        best_rel_B = []; 
        best_rel_id = 10;
        
        n_relaxations = size(Ps, 1); 
        set_min = zeros(n_relaxations, 1);  
        
        for i=1:n_relaxations
                P = Ps(i, :, :); 
                P = reshape(P, [5, 5]); 
                
                idx_col = discard_column(i); 
                
                set_indices = setdiff(2:7, idx_col+1);
                
                B = A(:, set_indices) * P;
                
                min_B = min(svd(B));
                set_min(i) = min_B; 
                
                if (min_B > best_min_sv)
                        best_min_sv = min_B; 
                        best_rel_B = B; 
                        best_rel_id = idx_col;
                end  %end update best relaxation
        end  % end for each relaxation

        struct_output = struct(); 
        struct_output.min_sv = best_min_sv;
        struct_output.best_rel = best_rel_B;
        struct_output.best_rel_id = best_rel_id;
        struct_output.set_min = set_min;
        
end  % end of function
