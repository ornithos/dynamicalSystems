function llh = calcLogLikelihood(obj)

    % Check for existence of Smoothed estimates
    if ~isfield(obj.posterior, 'smooth') || numel(obj.posterior.smooth.sigma) ~= obj.T
        fprintf('Smoothed estimates not found or incorrect format. Rerunning smoother...\n');
        obj = obj.posteriorFilter(false);
    end
    
    % get quantities
    P0          = obj.x0.sigma;
    m0          = obj.x0.mu;
    P           = obj.posterior.smooth.sigma;
    m           = obj.posterior.smooth.mu;
    G           = obj.posterior.smooth.G;
    y           = obj.y;
    T           = obj.T;
    
    Q           = obj.Q;
    R           = obj.R;
    H           = obj.H;
    A           = obj.A;
    
    [P0d, P0i]  = getDetAndInverse(P0);
    [Qd, Qi]    = getDetAndInverse(Q);
    [Rd, Ri]    = getDetAndInverse(R);
    
    SIGMA       = sum(cat(3,P{2:T}),3);
    m_Am1       = m(:,2:end) - A * m(:,1:end-1);
    C           = m_Am1 * m_Am1';
    y_Hm        = y - H * m;
    D           = y_Hm * y_Hm';
    F           = A * (SIGMA + P{1} - P{T}) * A';
    
    B     = 0;
    for tt = 2:T
        B = B + P{tt} * G{tt-1}';
    end
    B     = B*A';

    const       = zeros(2,1);
    const(1)    = (obj.d.x + obj.d.y) * log(2*pi);
    const(2)    = log(P0d);
    
    variable    = zeros(4 ,1);
    variable(1) = T*log(Qd) + T*log(Rd);
    variable(2) = trace(P0i * (P{1} + (m(:,1) - m0)*(m(:,1) - m0)'));
    variable(3) = trace(Qi  * (SIGMA - B - B' + F + C));
    variable(4) = trace(Ri  * (H * (SIGMA + P{1}) * H' + D));
    
    
    % ENTROPY??? WHERE ARE YOU ENTROPY?
    % -- this is not the true entropy, since the entropy is not marginally
    % independent between variables, but it gets us some of the way there.
    E     = T*(obj.d.x/2 + 0.5 * log(2*pi));
    for tt = 1:T
        E = E + log(det(P{tt}));
    end
    
    llh = - 0.5 * (sum(const) + sum(variable)); %+ E;
end

function [dd, ii] = getDetAndInverse(Z)
    [uu, ss, vv] = svd(Z);
    dd           = prod(diag(ss));
    ii           = uu * diag(1./(diag(ss))) * vv';
end