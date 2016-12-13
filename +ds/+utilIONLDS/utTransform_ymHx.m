function [ymHx, outerprod, XSP, Wc] = utTransform_ymHx(obj, alpha, beta, kappa)
    % out = utTransform_ymHx(obj, pars, alpha, beta, kappa)
    %
    % For a non-linear dynamical system object, get the Unscented Transform
    % of (y - h(x))
    
    assert(isa(obj, 'ds.dynamicalSystem'), 'first argument must be a valid dynamicalSystems object');
    n         = obj.d.x;
    d         = obj.d.y;
    emiParams = obj.par.emiNLParams;
    assert(isfield(emiParams, 'C'), 'parameter C not found in emission parameters');
    C         = emiParams.C;
    emiParams = emiParams.eta;
    assert(isnumeric(emiParams) && (all(size(emiParams) == [1, 4]) || all(size(emiParams) == [d, 4])), ...
        'emiParams.eta must be a matrix of size (1, 4) or (d.y, 4)');
    
    if nargin == 1
        alpha = 1;
        beta  = 0;
        kappa = 3-n;
    else
        assert(utils.is.numscal(alpha), 'alpha must be a numeric scalar');
        assert(utils.is.numscal(beta), 'beta must be a numeric scalar');
        assert(utils.is.numscal(kappa), 'kappa must be a numeric scalar');
    end
    
    %% get sigma point stuff
    
    lambda   = alpha^2 * (n + kappa) - n;
    scl      = sqrt(n + lambda);
    
    CX       = zeros(d, 2*n + 1, obj.d.T);
    XSP      = zeros(n, (2*n + 1)*obj.d.T);

    for tt = 1:obj.d.T
        sigma    = obj.infer.smooth.sigma{tt};
        mu       = obj.infer.smooth.mu(:,tt);
        
        % covariance spread
        S        = scl.*chol(sigma, 'lower');

        % sigma points
        spts     = [zeros(n,1), S, -S];
        spts     = bsxfun(@plus, spts, mu);
        
        % save
        CX(:,:,tt) = obj.par.emiNLParams.C * spts;
        XSP(:,(1+2*n)*(tt-1) + (1:2*n+1)) = spts;
    end
    
    % weights
    Wm       = ones(2*n+1, 1)./(2*(n + lambda));
    Wm(1)    = lambda/(n+lambda);
    Wc       = Wm;
    Wc(1)    = Wc(1) + (1 - alpha^2 + beta);
    
    %% Calculate outputs
    % run through function
    Hx       = ds.utilIONLDS.gen_sigmoid(CX, emiParams);
    ymHxSP   = bsxfun(@minus, permute(obj.y, [1, 3, 2]), Hx);  % tensor of Y - h(sigma points)

    ymHxSPWm = bsxfun(@times, ymHxSP, permute(Wm, [2,1,3]));   % weighted tensor ymHxSP (mean weight)
    
    if nargout <= 2
        ymHx     = squeeze(sum(ymHxSPWm,2));  % sum over weighted sigma points
    else
        ymHx     = reshape(ymHxSP, d, (1+2*n)*(obj.d.T));
    end

    % _______ Outer Product of (y - h(x))(y - h(x))' _____________________
    if nargout > 1
%         outerprod = zeros(d, d, obj.d.T);
        outerprod = zeros(d,d);
        ymHxSPWc = bsxfun(@times, ymHxSP, permute(Wc, [2,1,3]));   % weighted tensor ymHxSP (cov weight)
        
        for tt = 1:obj.d.T
%             outerprod(:,:,tt) = ymHxSPWc(:,:,tt) * ymHxSP(:,:,tt)';  % only need to weight one of the product
            outerprod = outerprod + ymHxSPWc(:,:,tt) * ymHxSP(:,:,tt)';  % only need to weight one of the product
        end
    end
end