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
    
    % get common interface for function whether they have parameters or
    % not.
    [f, Df, h, Dh] = obj.functionInterfaces;
    
    
    filterMu        = zeros(obj.d.x, obj.d.T);
    filterSigma     = cell(obj.d.T, 1);
    % initialise t-1 = 0 values to prior
    m               = obj.par.x0.mu;
    P               = obj.par.x0.sigma;
    d               = obj.d.y;
    Q               = obj.par.Q;
    R               = obj.par.R;
    
    % main forward step loop
    for tt = 1:obj.d.T
        F               = Df(m);
        m_minus         = f(m);
        P_minus         = F * P * F' + Q;

        H               = Dh(m_minus);
        S               = H * P_minus * H' + R;
        [Sinv, lam]     = utils.math.pinvAndEig(S, 1e-12);
        K               = (P_minus * H') * Sinv;
        deltaY          = obj.y(:,tt) - h(m_minus);
        m               = m_minus + K * deltaY;
        P               = P_minus - K * S * K';
        
        if bDoLLH
            llh         = llh - d*pi - 0.5*sum(log(lam)) - 0.5*deltaY'*Sinv*deltaY;
        end
        
        filterMu(:,tt)  = m;
        filterSigma{tt} = P;
    end
    
    obj.infer.filter.mu    = filterMu;
    obj.infer.filter.sigma = filterSigma;
    if bDoLLH; obj.infer.llh = llh; end;
    
    obj.infer.fpHash     = obj.parameterHash;
    obj.infer.fType      = 'EKF';
end