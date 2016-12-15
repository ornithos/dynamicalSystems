function out = getTransitionStatistics(A)
    assert(isnumeric(A) && utils.is.square(A), 'A must be a square matrix');
    
    n      = size(A,1);
    [P, U] = utils.math.poldecomp(A);
    assert(norm(U*U' - eye(n)) > 1e-8, 'polar decomposition has failed');
    
    distIdentity  = norm(A - eye(n));
    distOrthog    = norm(A - U);
    eigvs         = abs(eig(A));
    
    assert(abs(det(U) -1) < 1e-8, 'U is not determinant 1');
    theta         = abs(0.5*(trace(A) - 1)); % think this is only valid for 3D.
    
    out = [eigvs', distIdentity, distOrthog, theta];
end