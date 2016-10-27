function obj = smoothLinear(obj, bDoValidation)
   
    % obj = obj.smoothLinear
    % Perform posterior inference on each hidden state in the dynamical
    % system. This is the RTS Smoother updates / backwards equations for
    % lienar dynamics. Method of object type 'dynamicalSystem'. If the 
    % Kalman filter values do not exist, or exist for different parameter
    % values, we call filterKalman first.
    %
    % obj = obj.smoothLinear(false)
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
    if obj.infer.fpHash ~= obj.parameterHash
        fprintf('Filter not run or parameters changed. Rerunning filter...\n');
        obj = obj.posteriorFilter(false, false);
    end
    fMu          = obj.infer.filter.mu;
    fSigma       = obj.infer.filter.sigma;
    
    % initialise backward prior @ T = forward, and pre-allocate remaining.
    m            = fMu(:,obj.d.T);
    P            = fSigma{obj.d.T};
    smoothMu     = [zeros(obj.d.x, obj.d.T-1) m];
    smoothSigma  = vertcat(cell(obj.d.T-1, 1), P);
    smoothG      = cell(obj.d.T-1,1);
    
    A            = obj.par.A;
    Q            = obj.par.Q;
    H            = obj.par.H;
    R            = obj.par.R;
    
    % main forward step loop
    for tt = (obj.d.T-1):-1:1
        fP_t            = fSigma{tt};
        m_minus         = A * fMu(:,tt);
        P_minus         = A * fP_t * A' + Q;

%         barbS10         = obj.A * fP_t;
%         barbARev        = barbS10' * inv(P_minus);
%         barbSRev        = fP_t - barbARev * barbS10;
%         barb_m          = barbARev * m + fMu(:,tt) - barbARev * obj.A * fMu(:,tt);
%         barb_P          = barbARev * P * barbARev' + barbSRev;
        
        G               = (fP_t * A') / (P_minus);
        m               = fMu(:,tt) + G * (m - m_minus);
        P               = fP_t + G * (P - P_minus) * G';
        
        smoothMu(:,tt)  = m;
        smoothSigma{tt} = P;
        smoothG{tt} = G;
    end
    
    % x0 (purely to get G_0)
    fP_t            = obj.par.x0.sigma;
    m_minus         = A * obj.par.x0.mu;
    P_minus         = A * fP_t * A' + Q;
    G               = (fP_t * A') / (P_minus);
    m               = obj.par.x0.mu + G * (m - m_minus);
    P               = fP_t + G * (P - P_minus) * G';
    
    % save out
    obj.infer.smooth.mu    = smoothMu;
    obj.infer.smooth.sigma = smoothSigma;
    obj.infer.smooth.G     = smoothG;
    obj.infer.smooth.x0    = struct('mu', m, 'sigma', P, 'G', G);
    
    obj.infer.sType        = 'Linear';
end