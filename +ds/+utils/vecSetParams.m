function obj = vecSetParams(obj, v, forcePSD)
% obj = vecSetParams(obj, vec)
% take (LDS) vectorised parameters (in order A,Q,H,R) and set them in
% dynamicalSystems object, obj.

    assert(isa(obj, 'ds.dynamicalSystem'), 'argument must be a valid dynamicalSystems object');
    n = obj.d.x;
    d = obj.d.y;
    reqdEl = 2*n*n+n*d+d*d;
    assert(numel(v)==reqdEl, 'vec does not have the required number of parameters');
    
    obj.par.A            = reshape(v(1:n*n),n,n);
    Q                    = reshape(v(n*n+1:2*n*n),n,n);
    if forcePSD
        Q                = utils.math.projectPSD(Q, 1e-3);
    end 
    obj.par.Q            = Q;
    
    obj.par.H            = reshape(v(2*n*n+1:2*n*n+n*d),d,n);
    R                    = reshape(v(2*n*n+n*d+1:end),d,d);
    if forcePSD
        R                = utils.math.projectPSD(R, 1e-3);
    end
    obj.par.R            = R;
    
end