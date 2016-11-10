function obj = smooth(obj, sType, utpar, opts)
    % smooth(obj, bDoLLH, bDoValidation, utpar)
    % Dynamical System Smoothing for general dynamics as described in eg. 
    % Särkkä (2013). If the transition dynamics are non-linear, the
    % smoother will approximate the first two moments as required by 
    % either the first-order Taylor approximation (EKF), or the unscented
    % (cubature-type) transform as specified by fType. For smoothing, the
    % emission dynamics are not incorporated, so type of emission makes no
    % difference to the algorithm here. If the transition dynamics are 
    % linear, any fType specified will be ignored, and the relevant
    % dynamics calculated exactly as per the RTS equations
    %
    % INPUTS:
    % sType        - type of filter: 'Kalman', 'Extended', 'Unscented'.
    % utpar        - struct containing parameters for Unscented Transform
    %                (alpha, beta, kappa). Ignored if UKF not fType.
    % opts         - additional opts struct to modify internals. Includes:
    %                'bDoValidation': if false, switch off inp validation,
    %                'bIgnoreHash': neither check or perform paramater Hash
    %
    % OUTPUT:
    %  Returns object of type 'dynamicalSystem' with property
    %  infer.smooth as cell of means and variances of marginal posterior.
    %
    
    %% Argument defaults and validation
    
    if nargin < 4 || isempty(opts)
        opts = struct;
    end
    if nargin < 3 || isempty(utpar)
        utpar = struct;
    end
    if nargin < 2 || isempty(sType)
        if obj.evoLinear && obj.emiLinear
            sType = 'kal';
        else
            error(['smoother type (arg 1) not specified for nonlinear model. ', ...
                'Please specify algorithm.']);
        end
    end
    
    optsDefault  = struct('bDoValidation', true, 'bIgnoreHash', false);
    opts         = utils.struct.structCoalesce(opts, optsDefault);
    utparDefault = struct('alpha', 1e-3, 'beta', 2, 'kappa', 0);
    utpar        = utils.struct.structCoalesce(utpar, utparDefault);
    
    if opts.bDoValidation
        obj.validationInference;
    end
    
    bNumeric     = true;
    inpFtype     = sType;
    sType        = utils.filterTypeLookup(sType, bNumeric) - 1;
    if sType == 0 && ~(obj.evoLinear && obj.emiLinear)
        error('Unable to perform linear inference in non-linear model');
    elseif sType>0 && obj.evoLinear && obj.emiLinear
        warning(['Model is linear: exact inference will be performed using ', ...
                'RTS equations rather than %s type requested'], utils.filterTypeLookup(inpFtype));
    end
    
    % Check for existence of Filter
    if ~opts.bIgnoreHash && ~strcmp(obj.infer.fpHash, obj.parameterHash)
        fprintf('Filter not run or parameters changed. Rerunning filter...\n');
        obj = obj.filter(sType, false, utpar, opts);
    end
    
    %% Parameter setup
    fType1          = sType;
%     fType2          = fType;
    if obj.evoLinear; fType1 = 0; end
%     if obj.emiLinear; fType2 = 0; end
    
    parPredict   = getParams(obj, 1, fType1);
    fMu          = obj.infer.filter.mu;
    fSigma       = obj.infer.filter.sigma;
    
    % initialise backward prior @ T = forward, and pre-allocate remaining.
    m            = fMu(:,obj.d.T);
    P            = fSigma{obj.d.T};
    smoothMu     = [zeros(obj.d.x, obj.d.T-1) m];
    smoothSigma  = vertcat(cell(obj.d.T-1, 1), P);
    smoothG      = cell(obj.d.T-1,1);
    
    %% -- Main Loop --
    
    for tt = (obj.d.T-1):-1:1
        % Prediction step
        fP_tt           = fSigma{tt};
        fMu_tt          = fMu(:,tt);
        [m_minus, P_minus, covttp1] = assumedDensityTform(parPredict, fMu_tt, fP_tt, fType1, utpar);

        % Smoothing step
        G               = covttp1 / (P_minus);
        m               = fMu_tt + G * (m - m_minus);
        P               = fP_tt + G * (P - P_minus) * G';
        
        % save temporary
        smoothMu(:,tt)  = m;
        smoothSigma{tt} = P;
        smoothG{tt}     = G;
    end
    
    % x0 (purely to get G_0)
    fMu_tt          = obj.par.x0.mu;
    fP_tt           = obj.par.x0.sigma;
    [m_minus, P_minus, covttp1]  = assumedDensityTform(parPredict, fMu_tt, fP_tt, fType1, utpar);
   
    G               = covttp1 / (P_minus);
    m               = obj.par.x0.mu + G * (m - m_minus);
    P               = fP_tt + G * (P - P_minus) * G';
    
    % save out
    obj.infer.smooth.mu    = smoothMu;
    obj.infer.smooth.sigma = smoothSigma;
    obj.infer.smooth.G     = smoothG;
    obj.infer.smooth.x0    = struct('mu', m, 'sigma', P, 'G', G);
    
    obj.infer.sType        = 'UKF';
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