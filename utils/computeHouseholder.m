function R = computeHouseholder(x, e)
    x = x / norm(x, 2);
    e = e / norm(e, 2);
    
    u = x + sign(x(1)) * e;
    u = u / norm(u, 2);
    R = eye(3) - 2 * u * u';
end


