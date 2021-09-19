function struct_output = create2D2DCorrespondences(struct_input)

% CREATE2D2DCORRESPONDENCESMODEL: This function creates random correspondences
% A set of n_obs points are create inside a frustum (FOV AND min_depth, max_depth)
% we suppose the initial camera is at the origin 
% and create a random pose for the second camera 
% with bounded rotation (0.5 degrees) 
% and bounde translation (max_parallax) 
% The observations are computed as K [R, t] 
% where K is the intrinsic param. matrix 
% we do not include distortion 
% we add noise to the second set of observations wlog 
% and obtain the normalized vectors by K ^{-1}
% INPUT: 
%       - pt_number:       number of points
%       - noise:           noise in pixels for the observations 
%       - number_outliers: number of outliers you want to include
%       - FoV_par:         FOV in degrees
%       - max_parallax:    magnitude of the translation (in meters)
%       - min_depth, max_depth: depths (Z-axis) of the points 
%       - focal_length:    focal_length of the camera
%       - image_size:      image size of the camera as (rows, cols)

% OUTPUT: 
%       - obs1:            observations in the first frame 
%       - obs2:            observations in the second image 
%       - Ps:              3D world points coordiantes wrt the firts frame 
%       - t:               ground truth translation (with scale)
%       - R:               ground truth rotation 




% for tand fcn
FoV = mod(struct_input.FoV_par, 180); 



%% Generate random point-cloud
% inside the frustum defined by FOV and min/max depth
zpoints = struct_input.min_depth + (struct_input.max_depth - struct_input.min_depth) .* rand(1, struct_input.pt_number);
xmax = tand(FoV*0.5) * zpoints;
xmin = -tand(FoV*0.5) * zpoints;

xypoints = xmin + (xmax-xmin) .* rand(2, struct_input.pt_number);
points=[xypoints; zpoints];


%% generate random view-points

max_rotation = 0.5; % in degrees
position1 = zeros(3,1);  % frame1.t
rotation1 = eye(3);      % frame1.R

valid_position = false;
position2 = zeros(3,1); 
rotation2 = zeros(3, 1);
iter=0;

max_trials = 200;

while ((~valid_position) && (iter <= max_trials))
    rotation2 = generateBoundedR(max_rotation);

   
    % create random translation direction
    position2 = rand(3,1);  % always positive 
    % move the camera always backwards
    % scale
    position2 = struct_input.max_parallax .* position2 ./ norm(position2);
    
    rel_points = rotation2' * (points - position2);
    
    % all the points lie inside the FoV of the second camera ?
    
    % cehck FoV
    ratio_x = max(abs(rel_points(1) ./ rel_points(3)));  % max(tan(FoV/2))
    if ratio_x > tand(FoV*0.5)
        iter = iter + 1 ;
        % error: some points (X) were outside the FoV of the second camera
        continue
    end
    ratio_y = max(abs(rel_points(2) ./ rel_points(3)));  % max(tan(FoV/2))
    if ratio_y > tand(FoV*0.5)
        iter = iter + 1;
        % error: some points (Y) were outside the FoV of the second camera
        continue
    end
    
    % if you have arrived here, all the points are OK
    valid_position=true;  % not really needed
    break;
end % end- checking valid points & position


%% compute relative translation and rotation

R = rotation1' * rotation2;
t = rotation1' * (position2 - position1);


% intrinsic params
K = eye(3); 
K(1, 1) = struct_input.focal_length;
K(2, 2) = struct_input.focal_length;

% do not use this
% compute the size of the sensor with the FOV and focal length 
K(1, 3) = tand(FoV * 0.5) * struct_input.focal_length;
K(2, 3) = tand(FoV * 0.5) * struct_input.focal_length;



P0 = zeros(3, 4); P0(1:3, 1:3) = rotation1'; P0(1:3, end) = - rotation1' * position1;
P1 = zeros(3, 4); P1(1:3, 1:3) = rotation2'; P1(1:3, end) = - rotation2' * position2;

P0 = K * P0;
P1 = K * P1;


%% Now create the correspondences by looping through the cameras

obs1 = zeros(3, struct_input.pt_number);
obs2 = zeros(3, struct_input.pt_number);

for i=1:struct_input.pt_number
    % points wrt cameras
    body_point1 = points(:, i); 
    % body_point1 = rotation1' * (points(:,i)-position1);
    % body_point2 = rotation2' * (points(:,i)-position2);
    
    % homogeneous 
    body_point1 = [body_point1; 1];
    % body_point2 = [body_point2; 1];

    % project with pinhole camera model
    bearingVector1 = P0 * body_point1;
    bearingVector2 = P1 * body_point1;
    
    % homogenous 
    bearingVector1 = bearingVector1 ./ bearingVector1(end);
    bearingVector2 = bearingVector2 ./ bearingVector2(end);
    
    % add noise to the bearing vectors here
    bearingVector1_noisy = addNoiseGaussian(bearingVector1, struct_input.focal_length, struct_input.noise);
    bearingVector2_noisy = addNoiseGaussian(bearingVector2, struct_input.focal_length, struct_input.noise);

    % convert to mm
    mm_v1_noisy = inv(K) * bearingVector1_noisy;
    mm_v2_noisy = inv(K) * bearingVector2_noisy;
    
    obs1(:,i) = mm_v1_noisy ./norm(mm_v1_noisy, 2);
    obs2(:,i) = mm_v2_noisy ./norm(mm_v2_noisy, 2);
 
end
 
    
   
% save results 
% OUTPUT: 
%       - obs1:            observations in the first frame 
%       - obs2:            observations in the second image 
%       - Ps:              3D world points coordiantes wrt the firts frame 
%       - t:               ground truth translation (with scale)
%       - R:               ground truth rotation 

struct_output = Class2DObservationsOutput;
struct_output.obs1 = obs1; 
struct_output.obs2 = obs2;
struct_output.Ps = points;
struct_output.t = t;
struct_output.R = R;
struct_output.size_image = tand(FoV * 0.5) * struct_input.focal_length * 2;

