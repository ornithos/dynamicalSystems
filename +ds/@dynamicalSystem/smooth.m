function smooth(obj, sType, utpar, varargin)
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
    
    optsDefault  = struct('bDoValidation', true, 'bIgnoreHash', false, 'forceFilter', false, 'doLlh', false);
    opts         = utils.base.processVarargin(varargin, optsDefault);
    utparDefault = struct('alpha', 1e-3, 'beta', 2, 'kappa', 0);
    utpar        = utils.struct.structCoalesce(utpar, utparDefault);
    
    if opts.bDoValidation
        obj.validationInference;
    end
    
    bNumeric     = true;
    inpStype     = sType;
    sType        = ds.utils.filterTypeLookup(sType, bNumeric) - 1;
    if sType == 0 && ~(obj.evoLinear && obj.emiLinear)
        error('Unable to perform linear inference in non-linear model');
    elseif sType>0 && obj.evoLinear && obj.emiLinear
        warning(['Model is linear: exact inference will be performed using ', ...
                'RTS equations rather than %s type requested'], ds.utils.filterTypeLookup(inpStype));
    end
    
    % Check for existence of Filter
    if opts.forceFilter || opts.doLlh
        obj.filter(inpStype, opts.doLlh, utpar, opts);
    elseif ~opts.bIgnoreHash && obj.parametersChanged
        if obj.opts.verbose; fprintf('Filter not run or parameters changed. Rerunning filter...\n'); end
        obj.filter(inpStype, false, utpar, opts);
    end
    if ~strcmpi(obj.infer.fType, inpStype) && ~(strcmp(obj.infer.fType, 'Kalman') || ...
            sType==0)
        if obj.opts.warnings; warning('Different filter run to the specified smoother...\n'); end
    end
    
    %% Parameter setup
    fType1          = sType;
%     fType2          = fType;
    if obj.evoLinear; fType1 = 0; end
%     if obj.emiLinear; fType2 = 0; end
    
    parPredict   = ds.internal.getParams(obj, 1, fType1);
    fMu          = obj.infer.filter.mu;
    fSigma       = obj.infer.filter.sigma;
    
    % initialise backward prior @ T = forward, and pre-allocate remaining.
    m            = fMu(:,obj.d.T);
    P            = fSigma{obj.d.T};
    smoothMu     = [zeros(obj.d.x, obj.d.T-1) m];
    smoothSigma  = vertcat(cell(obj.d.T-1, 1), P);
    smoothG      = cell(obj.d.T-1,1);
    
    hasControl   = any(obj.hasControl);
    
    %% if nans --> end of sequence, we need to exclude backward step through
    % these: they should only use forward step
    bwdstart     = obj.d.T-1;   % usual case: no NaNs
%     while all(isnan(obj.y(:,bwdstart+1)))
%         
%         % ---- All of this machinery to get G -----------------------
%         if hasControl
%             u_t = obj.u(:,bwdstart+1);   %tt=T-1 => u_{T} because [N(Ax_t-1 + Bu_t, Q)]
%         else
%             u_t = [];
%         end
%         
%         fP_tt           = fSigma{bwdstart};
%         fMu_tt          = fMu(:,bwdstart);
%         [~, P_minus, covttp1] = ds.utils.assumedDensityTform(parPredict, fMu_tt, fP_tt, u_t, fType1, utpar);
%         
%         smoothG{bwdstart}      = covttp1 / P_minus; % in suffStats * by (now) filtered cov, so will be A * P.
%         % ------------------------------------------------------------
%         smoothMu(:, bwdstart)  = fMu_tt;          % filtered mean
%         smoothSigma{bwdstart}  = fP_tt;           % filtered covariance
%         bwdstart               = bwdstart - 1;
%     end
    
    %% -- Main Loop --
    for tt = bwdstart:-1:1
        % Prediction step
        if hasControl
            u_t = obj.u(:,tt+1);   %tt=T-1 => u_{T} because [N(Ax_t-1 + Bu_t, Q)]
        else
            u_t = [];
        end
        
        fP_tt           = fSigma{tt};
        fMu_tt          = fMu(:,tt);
        [m_minus, P_minus, covttp1] = ds.utils.assumedDensityTform(parPredict, fMu_tt, fP_tt, u_t, fType1, utpar);

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
    if hasControl
        u_t             = obj.u(:,1);
    end
    fMu_tt          = obj.par.x0.mu;
    fP_tt           = obj.par.x0.sigma;
    [m_minus, P_minus, covttp1]  = ds.utils.assumedDensityTform(parPredict, fMu_tt, fP_tt, u_t, fType1, utpar);
   
    G               = covttp1 / (P_minus);
    m               = obj.par.x0.mu + G * (m - m_minus);
    P               = fP_tt + G * (P - P_minus) * G';
    
    % save out
    obj.infer.smooth.mu    = smoothMu;
    obj.infer.smooth.sigma = smoothSigma;
    obj.infer.smooth.G     = smoothG;
    obj.infer.smooth.x0    = struct('mu', m, 'sigma', P, 'G', G);
    obj.infer.smooth.utpar = utpar;
    
    obj.infer.sType        = inpStype;
end
