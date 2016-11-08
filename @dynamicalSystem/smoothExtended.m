function obj = smoothExtended(obj, bDoValidation)
   
    % obj = obj.smoothExtended
    % Perform posterior inference on each hidden state in the dynamical
    % system. This is the Extended RTS Smoother backwards equations for
    % non-linear dynamics. Method of object type 'dynamicalSystem'. If the 
    % filter values do not exist, or exist for different parameter
    % values, we call filterExtended first.
    %
    % obj = obj.smoothExtended(false)
    % As above, but performs no input validation.
    %
    % OUTPUT:
    %  Returns object of type 'dynamicalSystem' with property
    %  smooth as matrix of means and cell of variances of marginal posterior.
    
    if nargin <2 || isempty(bDoValidation)
        bDoValidation = true;
    end
    
    if bDoValidation
        obj.validationInference;
    end
    
    % Check for existence of Filter
    if isempty(obj.infer.fpHash) || strcmp(obj.infer.fpHash, obj.parameterHash)
        fprintf('Filter not run or parameters changed. Rerunning filter...\n');
        obj = obj.filterExtended(false, false);
    end
    if ~strcmp(obj.infer.fType, 'EKF')
        fprintf('EKF not run (filter: %s). Running EKF instead...\n', obj.infer.fType);
        obj = obj.filterUnscented(false, false);
    end
    fMu          = obj.infer.filter.mu;
    fSigma       = obj.infer.filter.sigma;
    
    % initialise backward prior @ T = forward, and pre-allocate remaining.
    m            = fMu(:,obj.d.T);
    P            = fSigma{obj.d.T};
    smoothMu     = [zeros(obj.d.x, obj.d.T-1) m];
    smoothSigma  = vertcat(cell(obj.d.T-1, 1), P);
    smoothG      = cell(obj.d.T-1,1);
    
    Q            = obj.par.Q;
    [f, Df, ~, ~] = obj.functionInterfaces;
    
    % main forward step loop
    for tt = (obj.d.T-1):-1:1
        fP_t            = fSigma{tt};
        
        F               = Df(fMu(:,tt));
        m_minus         = f(fMu(:,tt));
        P_minus         = F * fP_t * F' + Q;
        
        G               = (fP_t * F') / (P_minus);
        m               = fMu(:,tt) + G * (m - m_minus);
        P               = fP_t + G * (P - P_minus) * G';
        
        smoothMu(:,tt)  = m;
        smoothSigma{tt} = P;
        smoothG{tt} = G;
    end
    
    % x0 (purely to get G_0)
    fP_t            = obj.par.x0.sigma;
    F               = Df(obj.par.x0.mu);
    m_minus         = f(obj.par.x0.mu);
    P_minus         = F * fP_t * F' + Q;
    G               = (fP_t * F') / (P_minus);
    m               = obj.par.x0.mu + G * (m - m_minus);
    P               = fP_t + G * (P - P_minus) * G';
    
    % save out
    obj.infer.smooth.mu    = smoothMu;
    obj.infer.smooth.sigma = smoothSigma;
    obj.infer.smooth.G     = smoothG;
    obj.infer.smooth.x0    = struct('mu', m, 'sigma', P, 'G', G);
    
    obj.infer.sType        = 'EKF';
end