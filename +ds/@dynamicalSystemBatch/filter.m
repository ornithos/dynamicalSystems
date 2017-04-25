function D = filter(obj, fType, bDoLLH, utpar, opts)
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
    %                update. In non-linear models this is the approximate
    %                assumed log likelihood.
    % utpar        - struct containing parameters for Unscented Transform
    %                (alpha, beta, kappa). Ignored if UKF not fType.
    % opts         - additional opts struct to modify internals. Includes:
    %                'bDoValidation': if false, switch off inp validation,
    %                'bIgnoreHash': neither check or perform paramater Hash
    %                'bCollectGradient' is currently implemented only as a
    %                trial for parameter A. It does not work with missing
    %                values as yet.
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
    
    optsDefault  = struct('bDoValidation', true, 'bIgnoreHash', false, ...
                    'bCollectGradient', false, 'T', obj.d.T);
    opts         = utils.struct.structCoalesce(opts, optsDefault);
    utparDefault = struct('alpha', 1, 'beta', 0, 'kappa', 0);
    utpar        = utils.struct.structCoalesce(utpar, utparDefault);
    
    if opts.bCollectGradient || nargout > 0
        error('gradient unavailable in batch mode at present');
    end
    
    if opts.bDoValidation
        obj.validationInference;
    end
    
    bNumeric     = true;
    inpFtype     = fType;
    fType        = ds.utils.filterTypeLookup(fType, bNumeric) - 1;
    if fType>0 && obj.evoLinear && obj.emiLinear
        warning(['Model is linear: exact inference will be performed using ', ...
                'Kalman equations rather than %s type requested'], ds.utils.filterTypeLookup(inpFtype));
    end
    
    %% Parameter set-up
    fType1          = fType;
    fType2          = fType;
    if obj.evoLinear; fType1 = 0; end
    if obj.emiLinear; fType2 = 0; end
    
    llh             = 0;
    d               = obj.d.y;
    if numel(opts.T) == 1; opts.T = repmat(opts.T, obj.d.n, 1); end
    T               = opts.T;

    fMuCell         = cell(obj.d.n,1);
    fSigmaCell      = cell(obj.d.n,1);
    
    fns             = obj.functionInterfaces(true);
    
    for nn = 1:obj.d.n
        % initialise t-1 = 0 values to prior
        
        parPredict      = ds.internal.getParams(obj, 1, fType1, fns, nn);
        parUpdate       = ds.internal.getParams(obj, 2, fType2, fns, nn);
    
        filterMu        = zeros(obj.d.x, T(nn));
        filterSigma     = cell(T(nn), 1);
        
        m               = obj.par.x0.mu{nn};
        P               = obj.par.x0.sigma{nn};
        
    %     u_t             = [];
    %     if any(obj.hasControl); u_t = zeros(obj.d.u,1); end

        %% -- Main Loop --

        for tt = 1:T(nn)

            % control signal is (t) for update and prediction in this version.
            if any(obj.hasControl)
                u_t = obj.u{nn}(:,tt);
            else
                u_t = [];
            end

            % Prediction step
            if obj.hasControl(1), u = u_t;
            else, u = [];
            end
            [m_minus, P_minus, ~] = ds.utils.assumedDensityTform(parPredict, m, P, u, fType1, utpar);
            P_minus               = (P_minus+P_minus')/2;

            % Update step -- deal with missing values here
            if obj.hasControl(2), u = u_t;
            else, u = [];
            end
            y                 = obj.y{nn}(:,tt);
            curMissing        = any(isnan(y));
            if ~curMissing
                [m_y, S, covxy]   = ds.utils.assumedDensityTform(parUpdate, m_minus, P_minus, u, fType2, utpar);
                S                 = (S+S')/2;
                K                 = covxy / S;
                %K                 = utils.math.mmInverseChol(covxy, cholS);
                deltaY            = y - m_y;
                m                 = m_minus + K*deltaY;
                P                 = P_minus - K * S * K';
            else
                [m, P, addlLLH]   = internal_updateStepMissing(obj, parUpdate, m_minus, P_minus, u, y, fType2, utpar, bDoLLH);
                if bDoLLH; llh = llh + addlLLH; end
            end

            if (bDoLLH || opts.bCollectGradient) && ~curMissing
                cholS       = cholcov(S);
                if bDoLLH
        %             [Sinv, lam] = utils.math.pinvAndEig(S, 1e-12);  % unstable if used in KF eqns (messed up!)
        %             llh         = llh - 0.5*d*log(2*pi) - 0.5*sum(log(lam)) - 0.5*deltaY'*Sinv*deltaY;
                    addlLLH     = utils.math.mvnlogpdfchol(deltaY', zeros(1,d), cholS);
        %             fprintf('diff in LLH calcs is: %.5e\n', addlLLH - log(mvnpdf(deltaY', zeros(1,d), S)));
                    llh         = llh + addlLLH;
                end
                if opts.bCollectGradient
                    D = ds.utils.dgradUpdate(obj, D, P_minus, S, K, deltaY, m, P, opts.bCollectGradient);
                end
            end
            % save
            filterMu(:,tt)  = m;
            filterSigma{tt} = P;
        end
        
        fMuCell{nn}    = filterMu;
        fSigmaCell{nn} = filterSigma;
    end
    
    obj.infer.filter.mu      = fMuCell;
    obj.infer.filter.sigma   = fSigmaCell;
    obj.infer.llh            = NaN;
    if bDoLLH; obj.infer.llh = llh; end
    if ~opts.bIgnoreHash; obj.infer.fpHash = obj.parameterHash; end
    obj.infer.fType          = inpFtype;
    obj.infer.sType          = '';
    obj.infer.filter.utpar   = utpar;
end

function [m, P, addlLLH] = internal_updateStepMissing(obj, parUpdate, m_minus, P_minus, u, y, fType2, utpar, bDoLlh)
    if all(isnan(y))
        m           = m_minus;
        P           = P_minus;
        addlLLH     = 0;
        return
    end
    
    % remove relevant indices
    missIdx       = isnan(y);
    Y             = y(~missIdx);
    
    if obj.emiLinear
        parUpdate.bias = parUpdate.bias(~missIdx);
        parUpdate.A    = parUpdate.A(~missIdx,:);            % actually H
        parUpdate.Q    = parUpdate.Q(~missIdx,~missIdx);     % actually R
    else
        [~,~,h,Dh]    = obj.functionInterfaces;
        % hack-y way to permit multiple statements in anonymous function
        AREF          = @(A,B) subsref(A, struct('type', '()', 'subs', {B}));
        parUpdate.f   = @(x) AREF(h(x), {find(~missIdx), ':'});
        parUpdate.Df  = @(x) AREF(Dh(x), {find(~missIdx), ':'});
        parUpdate.Q   = parUpdate.Q(~missIdx,~missIdx);     % actually R
    end
    
    % perform standard update
    [m_y, S, covxy]   = ds.utils.assumedDensityTform(parUpdate, m_minus, P_minus, u, fType2, utpar);
    S                 = (S+S')/2;
    K                 = covxy / S;
    %K                 = utils.math.mmInverseChol(covxy, cholS);
    deltaY            = Y - m_y;
    m                 = m_minus + K*deltaY;
    P                 = P_minus - K * S * K';

    if bDoLlh
        cholS       = cholcov(S);
        d           = sum(~missIdx);
        addlLLH     = utils.math.mvnlogpdfchol(deltaY', zeros(1,d), cholS);
    else
        addlLLH     = 0;
    end
end