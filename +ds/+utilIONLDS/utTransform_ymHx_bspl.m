function [ymHx, outerprod, XSP, Wc] = utTransform_ymHx_bspl(obj, alpha, beta, kappa)
    % out = utTransform_ymHx(obj, pars, alpha, beta, kappa)
    %
    % For a non-linear dynamical system object, get the Unscented Transform
    % of (y - h(x))
    
    assert(isa(obj, 'ds.dynamicalSystem'), 'first argument must be a valid dynamicalSystems object');
    n         = obj.d.x;
    d         = obj.d.y;
    emiParams = obj.par.emiNLParams;
    assert(isfield(emiParams, 'C'), 'parameter C not found in emission parameters');
    if ~isfield(emiParams, 'bias'); emiParams.bias = []; end
    numEtaReqd = numel(emiParams.bspl.t) - 2;
    assert(isnumeric(emiParams.eta) && all(size(emiParams.eta) == [d, numEtaReqd]), ...
        'emiParams.eta must be a matrix of size (d.y, %d)', numEtaReqd);
    assert(isnumeric(emiParams.bias) && (isempty(emiParams.bias) || all(size(emiParams.bias) == [d, 1])), ...
        'emiParams.bias must be a matrix of size (d.y, 1)');
    
    if nargin == 1
        alpha = 1;
        beta  = 0;
        kappa = 0;
    else
        assert(utils.is.scalarnum(alpha), 'alpha must be a numeric scalar');
        assert(utils.is.scalarnum(beta), 'beta must be a numeric scalar');
        assert(utils.is.scalarnum(kappa), 'kappa must be a numeric scalar');
    end
    
    %% get sigma point stuff
    anyNaN   = any(isnan(obj.y), 1);
    allNaN   = all(isnan(obj.y), 1);
    T        = find(~allNaN, 1, 'last');
    
    lambda   = alpha^2 * (n + kappa) - n;
    scl      = sqrt(n + lambda);
    
    CX       = zeros(d, 2*n + 1, T);
    Hx       = zeros(d, 2*n + 1, T);
    XSP      = zeros(n, (2*n + 1)*T);

    for tt = 1:T
        sigma    = obj.infer.smooth.sigma{tt};
        mu       = obj.infer.smooth.mu(:,tt);
        
        % covariance spread
        S        = scl.*chol(sigma, 'lower');

        % sigma points
        spts     = [zeros(n,1), S, -S];
        spts     = bsxfun(@plus, spts, mu);
        
        % save
        CX(:,:,tt) = obj.par.emiNLParams.C * spts;
        if ~isempty(emiParams.bias); CX(:,:,tt) = bsxfun(@plus, CX(:,:,tt), emiParams.bias); end
        Hx(:,:,tt) = ds.utilIONLDS.nonlinBsplNoC(CX(:,:,tt), emiParams.eta, emiParams.bspl);
        XSP(:,(1+2*n)*(tt-1) + (1:2*n+1)) = spts;
    end
    
    % weights
    Wm       = ones(2*n+1, 1)./(2*(n + lambda));
    Wm(1)    = lambda/(n+lambda);
    Wc       = Wm;
    Wc(1)    = Wc(1) + (1 - alpha^2 + beta);
    
    %% Calculate outputs
    % run through nonlinear function
    
    % spline does not work with tensors (but doesn't error either!!)
    % (moved into CX loop).
    
    ymHxSP   = bsxfun(@minus, permute(obj.y(:,1:T), [1, 3, 2]), Hx);  % tensor of Y - h(sigma points)

    ymHxSPWm = bsxfun(@times, ymHxSP, permute(Wm, [2,1,3]));   % weighted tensor ymHxSP (mean weight)
    
    if nargout <= 2
        ymHx     = squeeze(nansum(ymHxSPWm,2));  % sum over weighted sigma points
    else
        ymHx     = reshape(ymHxSP, d, (1+2*n)*T);
    end

    % _______ Outer Product of (y - h(x))(y - h(x))' _____________________
    if nargout > 1
        
%         outerprod = zeros(d, d, obj.d.T);
        outerprod = zeros(d,d);
        ymHxSPWc = bsxfun(@times, ymHxSP, permute(Wc, [2,1,3]));   % weighted tensor ymHxSP (cov weight)

        for tt = 1:T
%             if all(isnan(obj.y(:,tt))); continue; end
            curIncr   = ymHxSPWc(:,:,tt) * ymHxSP(:,:,tt)';  % weight (only) one of the product. 
            
            % handle missing values
            %    *---> cov(y_obs, y_unobs) = 0     (path blocked)
            %    *---> E(y_u - h_u(x_t))(y_u - h_u(x_t))') = R_u    since E(y_u) = h_u(x_t)
            if anyNaN(tt)
                if allNaN(tt)
                    curIncr = zeros(d,d);
                else
                    missing   = isnan(obj.y(:,tt));
                    curIncr(missing, ~missing) = 0;
                    curIncr(~missing, missing) = 0;
                    curIncr(missing,  missing) = obj.par.R(missing, missing);
                end
            end
            outerprod = outerprod + curIncr;
        end
    end
end