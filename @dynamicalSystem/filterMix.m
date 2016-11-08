function obj = filterMix(obj, fType, bDoValidation, utpar)
    % filterMix(obj, bDoLLH, bDoValidation, utpar)
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
    if nargin < 2 || isempty(fType)
        error('filter type (arg 2) not specified. Please specify algorithm.');
    end
    
    if bDoValidation
        obj.validationInference;
    end
    
    utparDefault = struct('alpha', 1, 'beta', 0, 'kappa', 0);
    utpar        = utils.struct.structCoalesce(utpar, utparDefault);
    
    switch lower(fType)
        case {'kalman', 'linear'}
            fType = 0;
        case {'extended','ekf'}
            fType = 1;
        case {'unscented', 'ukf', 'ut'}
            fType = 2;
        otherwise
            error('Unknown type given. Try ''linear'',''EKF'',or ''UKF''');
    end
    
    % get parameters
    fType1          = fType;
    fType2          = fType;
    if obj.evoLinear; fType1 = 0; end
    if obj.emiLinear; fType2 = 0; end
    
    parPredict      = getParams(obj, 1, fType1);
    parUpdate       = getParams(obj, 2, fType2);
    
    filterMu        = zeros(obj.d.x, obj.d.T);
    filterSigma     = cell(obj.d.T, 1);
    % initialise t-1 = 0 values to prior
    m               = obj.par.x0.mu;
    P               = obj.par.x0.sigma;
    
    %% -- Main Loop --
    
    for tt = 1:obj.d.T
        % Prediction step
        [m_minus, P_minus, ~] = assumedDensityTform(parPredict, m, P, fType1, utpar);

        % Filter step
        [m_y, S, covxy]   = assumedDensityTform(parUpdate, m_minus, P_minus, fType2, utpar);
        K                 = covxy / S;
        deltaY            = obj.y(:,tt) - m_y;
        m                 = m_minus + K*deltaY;
        P                 = P_minus - K * S * K';
        
        % save
        filterMu(:,tt)  = m;
        filterSigma{tt} = P;
    end
    
    
    obj.infer.filter.mu      = filterMu;
    obj.infer.filter.sigma   = filterSigma;
    obj.infer.llh            = NaN;
    obj.infer.fpHash         = obj.parameterHash;
    obj.infer.fType          = fType;
end

function [m, P, C] = assumedDensityTform(pars, m, P, type, utpar)
    % stage = 1: Prediction step
    % stage = 2: Update step

    % type = 0: Kalman (Linear)
    % type = 1: EKF
    % type = 2: UKF
    
    switch type
        case 0
            m         = pars.A * m;
            P         = pars.A * P * pars.A' + pars.Q;
            C         = P * pars.A';
        case 1
            F         = pars.Df(m);
            m         = pars.f(m);
            C         = P * F';
            P         = F * P * F' + pars.Q;
        case 2
            [m, P, C] = utils.unscentedTransform(pars.f, m, P, ...
                            utpar.alpha, utpar.beta, utpar.kappa);
            P         = P + pars.Q;
        otherwise
            error('Unknown assumed density type. Try 0=Kalman,1=EKF,2=UKF');
    end
end

function par = getParams(obj, stage, type)
    par = struct;
    [f,Df,h,Dh]    = obj.functionInterfaces;
    if stage == 1
            par.Q = obj.par.Q;
            if type == 0
                par.A  = obj.par.A;
            else
                par.f  = f;
                par.Df = Df;
            end
        elseif stage == 2
            par.Q = obj.par.R;
            if type == 0
                par.A = obj.par.H;
            else
                par.f  = h;
                par.Df = Dh;
            end
        else
            error('Unknown stage requested. Try 1=prediction or 2=update');
    end
end