function Rs = extractQuotientfromE(P1, P2, E)
    Rs = zeros(3, 6); 
    
    [R, t, ~, ~] = recoverPoseFromEssentialInliers(E, P1, P2);
    
    % compute householder
    R0=computeHouseholder(t, [0; 0; 1]);
    
    Rs(:, 1:3) = R0';
    Rs(:, 4:6) = R0*R;
end  % end -function
