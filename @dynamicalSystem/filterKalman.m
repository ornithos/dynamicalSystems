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
    m               = obj.par.x0.mu;
    P               = obj.par.x0.sigma;
    d               = obj.d.y;
    
    A               = obj.par.A;
    Q               = obj.par.Q;
    H               = obj.par.H;
    R               = obj.par.R;
    
    % main forward step loop
    for tt = 1:obj.d.T
        m_minus         = A * m;
        P_minus         = A * P * A' + Q;

        S               = H * P_minus * H' + R;
        [Sinv, lam]     = utils.math.pinvAndEig(S, 1e-12);
        K               = (P_minus * H') * Sinv;
        deltaY          = obj.y(:,tt) - H * m_minus;
        m               = m_minus + K * deltaY;
        P               = P_minus - K * S * K';
        
        if false && tt > obj.d.T - 5
            fprintf('%8f ', det(P));
        end
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
    obj.infer.fType      = 'Linear';
end