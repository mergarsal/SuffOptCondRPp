function [Rtrue, Ttrue, Rfalse, Tfalse] = recoverPoseFromEssential(E, P1, P2)
% RECOVERPOSEFROMESSENTIAL:
% This functions returns the rotation and translation terms
% given E and the set of correspondences points
% The cheirality constraint in forced without 
% actually computing the 3D Points
% 
% Method based on 
% Equivalent Constraints for Two-View Geometry: 
% Pose Solution/Pure Rotation Identification and 
% 3D Reconstruction
% Paper authors: Qi Cai, Yuanxin Wu, Lilian Zhang and Peike Zhang 



%% pick the physically realizable pose
[R1, R2, T1, T2]  = PoseEMat(E);

N=size(P1, 2);  % Here, N is the number of inliers since P1 is already filtered

% 1. Compute the right rotation X 
% 2. Compute the right translation Y
X = zeros(N, 9);
Y = zeros(N, 6);
for ii = 1:N
    p1 = P1(:, ii);
    p2 = P2(:, ii);
    X(ii, :) = kron(p2', p1');
    Y(ii, :) = [norm(p2, 2) * p1', norm(p1, 2)*p2'];
end
% 1. Compute the right rotation X 
ER1=E*E'*R1;
ER2=E*E'*R2;
M11 = X*ER1(:);
M12 = X*ER2(:);

if numel(find(M11 > 0)) >= numel(find(M12 > 0))
    Rtrue = R1;
    Rfalse = R2;
else
    Rtrue=R2;
    Rfalse=R1;
end

% 2. Compute the right translation Y
M21 = Y*[-Rtrue'; eye(3)]*T1;
M22 = Y*[-Rtrue'; eye(3)]*T2;
if numel(find(M21 > 0)) > numel(find(M22 > 0))
    Ttrue = T1;
    Tfalse = T2;
else
    Ttrue = T2;
    Tfalse = T1;
end


