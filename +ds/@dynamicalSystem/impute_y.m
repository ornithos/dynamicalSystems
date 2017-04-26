function [y, covY] = impute_y(obj, varargin)
    % y = obj.impute_y
    % impute missing values from data in dynamicalSystems object. This is
    % different from obj.getFittedValues in some sense only because fitted
    % values replace only the NaNs in y.
    %
    % y = obj.impute(y, 'filter', true)
    % y = obj.impute(y, 'smooth', true)
    % Choose either to use imputations only from filtered latent state
    % or smoothed (default).
    %
    % [y, covY] = obj.impute(...., 'variance', true)
    % Also return posterior covariance matrix over imputed Y. Typically
    % this will be singular (ie. positive *semi*-definite) because all
    % known coordinates will have zero (co)variance.
    
    y          = obj.y;
    missing    = isnan(y);
    anyMissing = any(missing, 1);
    
    if sum(sum(missing)) == 0
        if numel(varargin) == 0
            return
        end
    end
    
    % options
    optsDefault.filter    = false;     % neither assigned to ensure don't assign both with user input.
    optsDefault.smooth    = false;
    optsDefault.variance  = false;
    optsDefault.bIgnoreHash = false;   % do not check inference
    optsDefault.sample    = false;
    optsDefault.utpar     = struct('alpha', 1, 'beta', 0, 'kappa', 0);
    optsDefault.nRng      = 0;         % unused: for conformity with batch.
    
    opts                  = utils.base.processVarargin(varargin, optsDefault);
    if ~(opts.filter || opts.smooth)
        opts.smooth = true;             % default
    end
    if nargout > 1 && ~opts.variance
        warning('variance not specified in function call. Use ''variance'', ''true'' to retrieve variance');
        covY = [];
    end
    assert(xor(opts.filter, opts.smooth), 'Both filter and smooth are selected. Please select just one!');
    
    % retrieve filter / smooth mean 
    if opts.filter
        if ~opts.bIgnoreHash; obj.ensureInference('IMPUTE', 'filter'); end
        xhat = obj.infer.filter.mu;
        % (NOTE: UNLIKE in SUFFSTATS.m, THIS IS tt NOT tt + 1, AS NO X0)
        P    = obj.infer.filter.sigma;
    else
        if ~opts.bIgnoreHash; obj.ensureInference('IMPUTE', 'smooth'); end
        xhat = obj.infer.smooth.mu;
        % (NOTE: UNLIKE in SUFFSTATS.m, THIS IS tt NOT tt + 1, AS NO X0)
        P    = obj.infer.smooth.sigma;
    end
    
    % precompute cholesky if need samples
    cholR      = chol(obj.par.R);
    
    % UKF handles emissions differently as does not linearise: E(f(x)) ~= f(E(x))
    doUKF      = ~isempty(obj.infer.fType) && strcmpi(obj.infer.fType(1), 'u');
    
    % impute missing values
    for tt = find(anyMissing)
        u_t     = [];
        if obj.hasControl(2); u_t = obj.u(:,tt); end
        if ~opts.sample
            if ~doUKF
                yhat_tt  = obj.doEmission(xhat(:,tt), u_t);
            else
                yhat_tt  = obj.doEmission(xhat(:,tt), u_t, P{tt}, opts.utpar);
            end
        else
            tmpx     = xhat(:,tt) + chol(P{tt})' * randn(obj.d.x, 1);
            yhat_tt  = obj.doEmission(tmpx, u_t) + cholR' * randn(obj.d.y, 1);
        end
        mask     = missing(:,tt);
        y(mask,tt)  = yhat_tt(mask);
    end

    % get covariance
    if opts.variance
        d    = obj.d.y;
        
        if doUKF; [~,~,h,~]    = obj.functionInterfaces; end

        covY = repmat({zeros(d)}, obj.d.T,1);
        for tt = find(anyMissing)
            if ~doUKF
                covY_tt  = obj.par.H * P{tt} * obj.par.H' + obj.par.R;
            else
                utFns        = struct('f', h, 'Q', obj.par.R);
                u_t          = []; if obj.hasControl(2); u_t = obj.u(:,tt); end
                [~, covY_tt] = ds.utils.assumedDensityTform(utFns, xhat(:,tt), P{tt}, u_t, 2, opts.utpar);
            end
            mask     = missing(:,tt);
            covY_tt(~mask, ~mask) = 0;
            covY{tt} = covY_tt;
        end
    end 
end