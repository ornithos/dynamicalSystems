function obj = filter(obj, fType, bDoLLH, utpar, opts)
    % filter(obj, fType, bDoLLH, utpar, opts)
    % Dynamical System Filtering for general dynamics as described in eg. 
    % Särkkä (2013). If transition and/or emission dynamics are non-linear, 
    % the filter will approximate the first two moments as required by 
    % either the first-order Taylor approximation (EKF), or the unscented
    % (cubature-type) transform as specified by fType. If either/both the
    % dynamics are linear, any fType specified will be ignored, and the
    % relevant dynamics calculated exactly as per the Kalman filter
    % equations.
    %
    % INPUTS:
    % fType        - type of filter: 'Kalman', 'Extended', 'Unscented'.
    % bDoLLH       - Calculate the log likelihood alongside the filtering
    %                update. In non-linear models this is the approximate /
    %                assumed log likelihood.
    % utpar        - struct containing parameters for Unscented Transform
    %                (alpha, beta, kappa). Ignored if UKF not fType.
    % opts         - additional opts struct to modify internals. Includes:
    %                'bDoValidation': if false, switch off inp validation,
    %                'bIgnoreHash': neither check or perform paramater Hash
    %
    % OUTPUT:
    %  Returns object of type 'dynamicalSystem' with property
    %  infer.filter as cell of means and variances of marginal posterior.
    %
    
    %% Argument defaults and validation
    if nargin < 5 || isempty(opts)
        opts   = struct;
    end
    if nargin < 4 || isempty(utpar)
        utpar = struct;
    end
    if nargin < 3 || isempty(bDoLLH)
        bDoLLH = false;
    end
    if nargin < 2 || isempty(fType)
        if obj.evoLinear && obj.emiLinear
            fType = 'kal';
        else
            error(['filter type (arg 1) not specified for nonlinear model. ', ...
                'Please specify algorithm.']);
        end
    end
    
    optsDefault  = struct('bDoValidation', true, 'bIgnoreHash', false);
    opts         = utils.struct.structCoalesce(opts, optsDefault);
    utparDefault = struct('alpha', 1, 'beta', 0, 'kappa', 0);
    utpar        = utils.struct.structCoalesce(utpar, utparDefault);
    
    if opts.bDoValidation
        obj.validationInference;
    end
    
    bNumeric     = true;
    inpFtype     = fType;
    fType        = utils.filterTypeLookup(fType, bNumeric) - 1;
    if fType>0 && obj.evoLinear && obj.emiLinear
        warning(['Model is linear: exact inference will be performed using ', ...
                'Kalman equations rather than %s type requested'], utils.filterTypeLookup(inpFtype));
    end
    
    %% Parameter set-up
    fType1          = fType;
    fType2          = fType;
    if obj.evoLinear; fType1 = 0; end
    if obj.emiLinear; fType2 = 0; end
    
    parPredict      = getParams(obj, 1, fType1);
    parUpdate       = getParams(obj, 2, fType2);
    llh             = 0;
    d               = obj.d.y;
    
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
        [Sinv, lam]       = utils.math.pinvAndEig(S, 1e-12);
        K                 = covxy * Sinv;
        deltaY            = obj.y(:,tt) - m_y;
        m                 = m_minus + K*deltaY;
        P                 = P_minus - K * S * K';
        
        if bDoLLH
            llh         = llh - 0.5*d*log(2*pi) - 0.5*sum(log(lam)) - 0.5*deltaY'*Sinv*deltaY;
        end
        
        % save
        filterMu(:,tt)  = m;
        filterSigma{tt} = P;
    end
    
    
    obj.infer.filter.mu      = filterMu;
    obj.infer.filter.sigma   = filterSigma;
    obj.infer.llh            = NaN;
    if bDoLLH; obj.infer.llh = llh; end
    if ~opts.bIgnoreHash; obj.infer.fpHash = obj.parameterHash; end
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
            C         = P * pars.A';
            P         = pars.A * P * pars.A' + pars.Q;
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