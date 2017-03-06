function D = dgradInitialise(obj)
    assert(isa(obj, 'ds.dynamicalSystem'), 'input is not a dynamical system');
    
    m         = obj.d.y;
    n         = obj.d.x;
    
    %I_m       = sparse(1:m,1:m,1,m,m);
    I_n       = sparse(1:n,1:n,1,n,n);
    D         = struct;
    
    % gradient and summations
    D.G_A     = zeros(1,n^2);   % in order to facilitate k-step
    D.phi     = zeros(1,n^2);   % phi = G_A (e.g.) for 1-step
    
    % initialise m0, P0
    D.m_minus = kron(obj.par.x0.mu', I_n);
    D.P_minus = kron(obj.par.A * obj.par.x0.sigma, I_n);
    D.P_minus = utils.math.vecPermTranspose(n,n) * D.P_minus + D.P_minus;
    
    % pre-allocate
    D.v       = zeros(m, n^2);
    D.S       = zeros(n^2, n^2);
    D.K       = zeros(m*n, n^2);
    D.m       = zeros(n, n^2);
    D.P       = zeros(n^2, n^2);
end