function obj = posteriorSmooth(obj, bDoValidation)
   
    % obj = obj.posteriorSmooth
    % Perform posterior inference on each hidden state in the dynamical
    % system. This is the RTS Smoother updates / backwards equations.
    % Method of object type 'dynamicalSystem'. If the Kalman Smoother
    % updates do not exist, we call posteriorFilter first.
    %
    % obj = obj.posteriorSmooth(false)
    % As above, but performs no input validation.
    %
    % OUTPUT:
    %  Returns object of type 'dynamicalSystem' with property
    %  posterior.smooth as matrix of means and cell of variances of
    %  marginal posterior.
    
    if nargin <2 || isempty(bDoValidation)
        bDoValidation = true;
    end
    
    if bDoValidation
        obj.validationInference;
    end
    
    % Check for existence of Filter
    if ~isfield(obj.posterior, 'filter') || numel(obj.posterior.filter.sigma) ~= obj.T
        fprintf('Filter estimates not found or incorrect format. Rerunning filter...\n');
        obj = obj.posteriorFilter(false, false);
    end
    fMu          = obj.posterior.filter.mu;
    fSigma       = obj.posterior.filter.sigma;
    
    % initialise backward prior @ T = forward, and pre-allocate remaining.
    m            = fMu(:,obj.T);
    P            = fSigma{obj.T};
    smoothMu     = [zeros(obj.d.x, obj.T-1) m];
    smoothSigma  = vertcat(cell(obj.T-1, 1), P);
    smoothG      = cell(obj.T-1,1);
    
    
    % main forward step loop
    for tt = (obj.T-1):-1:1
        fP_t            = fSigma{tt};
        m_minus         = obj.A * fMu(:,tt);
        P_minus         = obj.A * fP_t * obj.A' + obj.Q;

%         barbS10         = obj.A * fP_t;
%         barbARev        = barbS10' * inv(P_minus);
%         barbSRev        = fP_t - barbARev * barbS10;
%         barb_m          = barbARev * m + fMu(:,tt) - barbARev * obj.A * fMu(:,tt);
%         barb_P          = barbARev * P * barbARev' + barbSRev;
        
        G               = (fP_t * obj.A') / (P_minus);
        m               = fMu(:,tt) + G * (m - m_minus);
        P               = fP_t + G * (P - P_minus) * G';
        
        smoothMu(:,tt)  = m;
        smoothSigma{tt} = P;
        smoothG{tt} = G;
    end
    
    % x0 (purely to get G_0)
    fP_t            = obj.x0.sigma;
    m_minus         = obj.A * obj.x0.mu;
    P_minus         = obj.A * fP_t * obj.A' + obj.Q;
    G               = (fP_t * obj.A') / (P_minus);
    m               = obj.x0.mu + G * (m - m_minus);
    P               = fP_t + G * (P - P_minus) * G';
    
    % save out
    obj.posterior.smooth.mu    = smoothMu;
    obj.posterior.smooth.sigma = smoothSigma;
    obj.posterior.smooth.G     = smoothG;
    obj.posterior.smooth.x0    = struct('mu', m, 'sigma', P, 'G', G);
end