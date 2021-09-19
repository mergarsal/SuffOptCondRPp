function E = extractEFromQuotient(X)
    % from manopt: 
    % function: essential_costE2cost
    e3hat = [0 -1 0; 1 0 0; 0 0 0];

    RA = X(:,1:3,:); 
    RB = X(:,4:6,:); 
    E = multiprod(multiprod(multitransp(RA), e3hat), RB); 
end

