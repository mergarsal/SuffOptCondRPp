clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    Example sufficient optimality condition   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('shuffle');


%% Parameters for problem generation
noise = 0.5; 
nr_correspondences = 100; 
max_parallax = 2.0; 
min_depth=1.0;
max_depth=8.0;
focal_l = 800; 
Weights = ones(nr_correspondences, 1); 
fov = 100;


% create matrices P
Ps = createScaleMatrices();
% relaxation indices
array_idx_relaxations = 1:6; 
THRESHOLD_COND = 2e-04;
    
    
% 1. Generate data
struct_input = Class2DObservationsInput(); 
% define the params 
struct_input.pt_number = nr_correspondences;
struct_input.noise = noise;          % in pixels 
struct_input.FoV_par = fov;                % in degrees
struct_input.max_parallax = max_parallax;     % in meters
struct_input.min_depth = min_depth;           % in meters
struct_input.max_depth = max_depth;           % in meters                  
struct_input.focal_length = focal_l;             % in pixels, if 0 => compute from image size 
% Create problem 
struct_output = create2D2DCorrespondences(struct_input);

% extract data 
P1 = struct_output.obs1;  % observations in camera 0
P2 = struct_output.obs2;  % observations in camera 1
Rgt = struct_output.R;    % relative rotation
tgt = struct_output.t;    % relative translation
PPs = struct_output.Ps;   % 3D world points
tgt = tgt ./ norm(tgt);
% Compute the essential matrix as  E = [t]_x R
Egt = [[0 -tgt(3) tgt(2)]; [tgt(3) 0 -tgt(1)]; [-tgt(2) tgt(1) 0]]*Rgt;


% Create data matrix
C = constructCoeffMatrixC(P1, P2,  Weights);
% padded with zeros
Q = zeros(12); Q(1:9, 1:9) = C;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%       SMALLEST EIGENVALUE        %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the eigenvector associated with 
% the smallest eigenvalue of the data matrix
[uc, ~]=eig(C);
ecc = uc(:, 1);
ecc = ecc .* sqrt(2);
% estimate its left singular vector 
% with smallest singular value
Ecc = reshape(ecc, [3, 3]);
[ucc, dcc, vcc]= svd(Ecc);
tcc = ucc(:, end); 
xcc = [ecc; tcc];


% Estimate the best relaxation for this problem
sigma_struct = computeMinSingularValue(xcc, Ps, array_idx_relaxations);
% Get the data for the best relaxation
sigma_cc = sigma_struct.min_sv;
best_relaxation = sigma_struct.best_rel_id; 
% select the matrix P associated to the best relaxation
P = Ps(best_relaxation, :, :); 
P = reshape(P, [5, 5]); 

    
%%  PARAMETERS FOR THIS PROBLEM INSTANCE  %%
% create constraint matrices, including A1 = blkdiag(zeros(9), eye(3));
AsR_all = createAllConstraintMatricesReduced();
% compute the indices of the matrices
% we remove the best relaxation and the first matrix
idx_relaxations = setdiff(2:7, best_relaxation+1);
% retrieve constraint matrices for this relaxation
AsR = AsR_all(:, :, idx_relaxations); 
% pre-prend the matrix A1
AsR_complete = AsR_all(:, :, [1, idx_relaxations]); 

% create constraint wide matrix 
As = reshape(AsR, [12, 12*5]);
As_complete = [AsR_all(:, :, 1), As]; 



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Estimation solution
% we use the 8-pt algorithm as initial guess
E_initial = ucc * diag([1, 1, 0]) .* sqrt(2) *vcc'; 
% convert the initial matrix to the param. employed by the solver
Rs = extractQuotientfromE(P1, P2, E_initial);

% 3.2. Refine initial guess with iterative method (optimization on manifold)
out_mani = solveRPpManifold(C, Rs);
% retrieve results
x = out_mani.x; 
E_mani = out_mani.E; 
f_mani = out_mani.f;   % cost of the solution


%% check optimality with certifier first
% compute translation associated with E
e_mani = E_mani(:);
[ue, de, ve] = svd(E_mani);
t_mani = ue(:, 3);
xmani = [E_mani(:); t_mani]; 

% Parameters for the certifier
in_cert = struct(); 
in_cert.e_min_svd = THRESHOLD_COND;
in_cert.Q = Q; 
in_cert.E = E_mani; 
in_cert.t = t_mani; 
in_cert.As = AsR;
in_cert.f_hat  = f_mani;
% Call certifier
out_cert = checkOptimalityCertifier(in_cert);


%% Sufficient condition
struct_input_cond_mani = struct(); 
% INPUT: 
% - `threshold`:          threshold epsilon in paper
% - `min_sv_relaxation`:  min. singular value of the chosen relaxation
% - `xcc':                eigenvector associated with the smallest eigenvalue of Q
% - `Bs`:                 Wide matrix with all the constraint matrices for this relaxation: 12 x (12 ?? 5)
% - `Q`:                  data matrix Q (12 x 12)
% - `x`:                  solution to the relaxtive pose problem as [e, t]
% - `f`:                  objective value attain by the solution x

struct_input_cond_mani.threshold = THRESHOLD_COND;  
struct_input_cond_mani.min_sv_relaxation = sigma_cc;
struct_input_cond_mani.xcc = xcc;
struct_input_cond_mani.Bs = As;
struct_input_cond_mani.Q = Q;
struct_input_cond_mani.P = P;


%% Check solution from on-manifold optimization with condition
struct_input_cond_mani.x = xmani;
struct_input_cond_mani.f = f_mani;
% OUTPUT
% - `approx_mu_min`:      approx. of the minimum eigenvalue of M 
% - `is_opt`:             flag: optimality could be certified by bounds

struct_output_cond_mani = checkOptimalityCondition(struct_input_cond_mani);




