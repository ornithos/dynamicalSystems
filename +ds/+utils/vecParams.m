function v = vecParams(obj)
% vecParams(obj)
%  extract parameters from dynamicalSystem object and 'vec' them (stack
%  them all in a vector.
%  This is currently only available for linear dynamicalSystems.

    assert(isa(obj, 'ds.dynamicalSystem'), 'argument must be a valid dynamicalSystems object');
    assert(obj.emiLinear && obj.evoLinear, 'dynamicalSystems object must be linear!');
    n = obj.d.x;
    d = obj.d.y;
    
    v                    = NaN(2*n*n+n*d+d*d,1);
    v(1:n*n)             = obj.par.A(:);
    v(n*n+1:2*n*n)       = obj.par.Q(:);
    v(2*n*n+1:2*n*n+n*d) = obj.par.H(:);
    v(2*n*n+n*d+1:end)   = obj.par.R(:);
end