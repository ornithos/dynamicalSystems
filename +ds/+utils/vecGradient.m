function v = vecGradient(obj)
% vecGradient(obj)
%  calculate the gradient of a dynamicalSystems object at the current
%  parameter values, and 'vec' the result (stack the result in a vector).
%  (currently only coded for linear case).

    assert(isa(obj, 'ds.dynamicalSystem'), 'argument must be a valid dynamicalSystems object');
    n = obj.d.x;
    d = obj.d.y;
    G = obj.getGradient;
    
    v                    = NaN(2*n*n+n*d+d*d,1);
    v(1:n*n)             = G.A(:);
    v(n*n+1:2*n*n)       = G.Q(:);
    v(2*n*n+1:2*n*n+n*d) = G.H(:);
    v(2*n*n+n*d+1:end)   = G.R(:);
end