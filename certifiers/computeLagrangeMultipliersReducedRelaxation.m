function [L_hat, A_x, b, res]=computeLagrangeMultipliersReducedRelaxation(Q, As, e_opt, t)
%(Q-L*) * x* = 0

x_opt=[e_opt; t];
b = Q*x_opt;
x1=x_opt(1);x2=x_opt(2);x3=x_opt(3);x4=x_opt(4);x5=x_opt(5);x6=x_opt(6);
x7=x_opt(7);x8=x_opt(8);x9=x_opt(9);x10=x_opt(10);x11=x_opt(11);x12=x_opt(12);

% create matrix A wint coeff.
[n_r, n_c, n_mult] = size(As); 

As_stack = reshape(As, [12, 12*n_mult]); 

X = kron(eye(n_mult), x_opt); 
A_x = As_stack * X;

[Q, R]=qr(A_x);
% economy decomposition of A
% [Q, R]=qr(A, 0);
% Solve A*L_hat = b via QR fact.
L_hat = pinv(R)*(Q'*b);


% eig(A'*A);
res = norm(A_x*L_hat-b);

% least-square
[L_hat, flag, resrel, iter, resvec] = lsqr(A_x, b, []);


