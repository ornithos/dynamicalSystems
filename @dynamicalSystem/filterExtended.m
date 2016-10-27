function obj = filterExtended(obj, bDoLLH, bDoValidation)
   
    % obj = obj.filterExtended
    % Perform posterior inference on each hidden state in the dynamical
    % system. This is the Extended Kalman Filter (EKF) forward eqns for
    % non-linear dynamics. Method of object type 'dynamicalSystem'.
    %
    % obj = obj.filterExtended(true)
    % As above, but returns Log Likelihood.
    %
    % obj = obj.filterExtended([], false)
    % As above, but performs no input validation.
    %
    % OUTPUT:
    %  Returns object of type 'dynamicalSystem' with property
    %  filter as cell of means and variances of marginal posterior.
    
    % FOR ALEX:
    % * How does LLH work here - can we calculate it in the same way?
    % * How do we record which type of inference was used?
    
    if nargin <3 || isempty(bDoValidation)
        bDoValidation = true;
    end
    
    if nargin < 2 || isempty(bDoLLH)
        bDoLLH = false;
    end
    if bDoLLH
        llh = 0;
    end
    
    if bDoValidation
        obj.validationInference;
    end
    
    filterMu        = zeros(obj.d.x, obj.d.T);
    filterSigma     = cell(obj.d.T, 1);
    % initialise t-1 = 0 values to prior
    m               = obj.x0.mu;
    P               = obj.x0.sigma;
    d               = obj.d.y;
    
    % main forward step loop
    for tt = 1:obj.d.T
        F               = obj.Df(m, obj.evoNLParams);
        m_minus         = obj.f(m, obj.evoNLParams);
        P_minus         = F * P * F' + obj.Q;

        H               = obj.Dh(m_minus, obj.emiNLParams);
        S               = H * P_minus * H' + obj.R;
        [Sinv, lam]     = utils.math.pinvAndEig(S, 1e-12);
        K               = (P_minus * H') * Sinv;
        m               = m_minus + K * (obj.y(:,tt) - obj.h(m_minus, obj.emiNLParams));
        P               = P_minus - K * S * K';
        
        if bDoLLH
            deltaY      = obj.y(:,tt) - obj.h(m_minus, obj.emiNLParams);
            llh         = llh - d*pi - 0.5*sum(log(lam)) - 0.5*deltaY'*Sinv*deltaY;
        end
        
        filterMu(:,tt)  = m;
        filterSigma{tt} = P;
    end
    
    obj.filter.mu    = filterMu;
    obj.filter.sigma = filterSigma;
    if bDoLLH; obj.llh = llh; end;
    
    obj.fpHash     = obj.parameterHash;
end