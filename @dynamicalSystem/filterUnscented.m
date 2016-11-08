function obj = filterUnscented(obj, bDoLLH, bDoValidation, utpar)
    % filterUnscented(obj, bDoLLH, bDoValidation, utpar)
    % Unscented Kalman Filter for non-linear dynamics as described in eg. 
    % Wan and Van Der Merwe (2000). This is algorithmically equivalent to
    % the Kalman filter except the expectations and covariances are
    % calculated using the Unscented Transform instead of exactly in the
    % linear case.
    %
    % INPUTS:
    % fn     - a function handle capable of ingesting a matrix (where each
    %          datapoint corresponds to a column) and outputting a matrix
    %          (where each column corresponds to the mapped datapoint).
    % mu     - The mean of the initial Gaussian distribution. Must be of
    %          same dimension as number of rows in the matrix described
    %          above.
    % sigma  - The covariance matrix of the Gaussian distribution.
    %
    % OPTIONAL:
    % All the below parameters may be specified or none of them.
    % alpha  - scaling parameter (related to how clustered the sigma points
    %          are around the mean.) Suggest 10^-3
    % beta   - adjustment made for kurtosis of transformed distn. Suggest 2
    % lambda - original scaling parameter of UT. Suggest (3 - dim(mu))
    %
    % OUTPUTS:
    % mu     -  Mean of the closest transformed distribution according to
    %          the Unscented transform.
    % sigma  -  Covariance matrix of the same distribution.
    %
    
    %% Argument defaults and validation
    if nargin < 4 || isempty(utpar)
        utpar = struct;
    end
    
    if nargin <3 || isempty(bDoValidation)
        bDoValidation = true;
    end
    
    if nargin < 2 || isempty(bDoLLH)
        bDoLLH = false;
    end
    llh = 0;
    
    if bDoValidation
        obj.validationInference;
    end
    
    [f,~,h,~]    = obj.functionInterfaces;
    utparDefault = struct('alpha', 1, 'beta', 0, 'kappa', 0);
    utpar        = utils.struct.structCoalesce(utpar, utparDefault);
    
    
    
    filterMu        = zeros(obj.d.x, obj.d.T);
    filterSigma     = cell(obj.d.T, 1);
    % initialise t-1 = 0 values to prior
    m               = obj.par.x0.mu;
    P               = obj.par.x0.sigma;
    Q               = obj.par.Q;
    R               = obj.par.R;
    
    %% -- Main Loop --
    
    for tt = 1:obj.d.T
        % Prediction step
        [m_minus, c_m] = utils.unscentedTransform(f, m, P, utpar.alpha, ...
                            utpar.beta, utpar.kappa);
        P_minus        = c_m + Q;

        % Filter step
        [m_y, c_y, covxy] = utils.unscentedTransform(h, m_minus, P_minus, ...
                            utpar.alpha, utpar.beta, utpar.kappa);
        S                 = c_y + R;
        [Sinv, lam]       = utils.math.pinvAndEig(S, 1e-12);
        K                 = covxy * Sinv;
        deltaY            = obj.y(:,tt) - m_y;
        m                 = m_minus + K*deltaY;
        P                 = P_minus - K * S * K';
    
        if bDoLLH
            llh         = llh - d*pi - 0.5*sum(log(lam)) - 0.5*deltaY'*Sinv*deltaY;
        end
        
        filterMu(:,tt)  = m;
        filterSigma{tt} = P;
    end
    
    
    obj.infer.filter.mu      = filterMu;
    obj.infer.filter.sigma   = filterSigma;
    if bDoLLH; obj.infer.llh = llh; end;
    
    obj.infer.fpHash         = obj.parameterHash;
    obj.infer.fType          = 'UKF';
end