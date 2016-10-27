function obj = filterKalman(obj, bDoLLH, bDoValidation)
   
    % obj = obj.filterKalman
    % Perform posterior inference on each hidden state in the dynamical
    % system. This is the Kalman Filter updates / forward equations for
    % linear dynamics. Method of object type 'dynamicalSystem'.
    %
    % obj = obj.filterKalman(true)
    % As above, but returns Log Likelihood.
    %
    % obj = obj.filterKalman([], false)
    % As above, but performs no input validation.
    %
    % OUTPUT:
    %  Returns object of type 'dynamicalSystem' with property
    %  filter as cell of means and variances of marginal posterior.
    
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
        m_minus         = obj.A * m;
        P_minus         = obj.A * P * obj.A' + obj.Q;

        S               = obj.H * P_minus * obj.H' + obj.R;
        [Sinv, lam]     = utils.math.pinvAndEig(S, 1e-12);
        K               = (P_minus * obj.H') * Sinv;
        m               = m_minus + K * (obj.y(:,tt) - obj.H * m_minus);
        P               = P_minus - K * S * K';
        
        if false && tt > obj.d.T - 5
            fprintf('%8f ', det(P));
        end
        if bDoLLH
            deltaY      = obj.y(:,tt) - obj.H * m_minus;
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