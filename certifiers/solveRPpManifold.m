function str_out = solveRPpManifold(Q, Rs)


% 2. Cretae manifold
manifold = essentialfactory(1, 'unsigned');
problem.M = manifold;

% 3. Define cost
problem.cost = @(X) cost(X, Q);
costE = @(E) trace(vec(E)'*Q*vec(E));

problem.cost = @cost;
problem.egrad = @egrad;
problem.ehess = @ehess;

egradE = @(E) reshape(2 * Q * E(:), [3,3]);

ehessE = @(E, U) reshape(2 * Q * U(:), [3, 3]);


    function f = cost(X)
       f = essential_costE2cost(X, costE); 
    end
   
    function g = egrad(X)
        g = essential_egradE2egrad(X, egradE); % Converts gradient in E to X.
    end

    function gdot = ehess(X, S)
        gdot = essential_ehessE2ehess(X, egradE, ehessE, S); % Converts Hessian in E to X.
    end
    

 

options.tolgradnorm = 1e-10;  
Hes=hessianmatrix(problem, Rs);
options.maxiters = 100;

[x, xcost, info, options] = trustregions(problem, Rs, options);

E = extractEFromQuotient(x);

str_out = struct(); 
str_out.x = x; 
str_out.E = E; 
str_out.f = xcost; 
str_out.info = info; 
str_out.options = options;
end

