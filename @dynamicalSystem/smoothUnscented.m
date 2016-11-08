function obj = smoothUnscented(obj, bDoValidation, utpar)
    % obj = obj.smoothUnscented
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
    
    if nargin < 3 || isempty(utpar)
        utpar = struct;
    end
    
    if nargin <2 || isempty(bDoValidation)
        bDoValidation = true;
    end
    
    if bDoValidation
        obj.validationInference;
    end
    
    % Check for existence of Filter
    if ~strcmp(obj.infer.fpHash, obj.parameterHash)
        fprintf('Filter not run or parameters changed. Rerunning filter...\n');
        obj = obj.filterUnscented(false, false);
    end
    
    if ~strcmp(obj.infer.fType, 'UKF')
        fprintf('UKF not run (filter: %s). Running UKF instead...\n', obj.infer.fType);
        obj = obj.filterUnscented(false, false);
    end
    
    utparDefault = struct('alpha', 1e-3, 'beta', 2, 'kappa', 0);
    utpar        = utils.struct.structCoalesce(utpar, utparDefault);
    
    fMu          = obj.infer.filter.mu;
    fSigma       = obj.infer.filter.sigma;
    
    % initialise backward prior @ T = forward, and pre-allocate remaining.
    m            = fMu(:,obj.d.T);
    P            = fSigma{obj.d.T};
    smoothMu     = [zeros(obj.d.x, obj.d.T-1) m];
    smoothSigma  = vertcat(cell(obj.d.T-1, 1), P);
    smoothG      = cell(obj.d.T-1,1);
    
    Q            = obj.par.Q;
    [f,~,~,~]    = obj.functionInterfaces;
    
    %% -- Main Loop --
    
    for tt = (obj.d.T-1):-1:1
        fP_t            = fSigma{tt};
        fMu_t           = fMu(:,tt);
        [m_minus, c_m, covttp1]  = utils.unscentedTransform(f, fMu_t, fP_t,  ...
                                    utpar.alpha, utpar.beta, utpar.kappa);
        P_minus         = c_m + Q;
        
        G               = covttp1 / (P_minus);
        m               = fMu(:,tt) + G * (m - m_minus);
        P               = fP_t + G * (P - P_minus) * G';
        
        smoothMu(:,tt)  = m;
        smoothSigma{tt} = P;
        smoothG{tt}     = G;
    end
    
    % x0 (purely to get G_0)
    fP_t            = obj.par.x0.sigma;
    [m_minus, c_m, covttp1]  = utils.unscentedTransform(f, m, fP_t, utpar.alpha, ...
                                    utpar.beta, utpar.kappa);
    P_minus         = c_m + Q;
    G               = covttp1 / (P_minus);
    m               = obj.par.x0.mu + G * (m - m_minus);
    P               = fP_t + G * (P - P_minus) * G';
    
    % save out
    obj.infer.smooth.mu    = smoothMu;
    obj.infer.smooth.sigma = smoothSigma;
    obj.infer.smooth.G     = smoothG;
    obj.infer.smooth.x0    = struct('mu', m, 'sigma', P, 'G', G);
    
    obj.infer.sType        = 'UKF';
end