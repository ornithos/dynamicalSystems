function [f, grad, more] = bsplineGrad(obj, x, utpar, varargin)
    % admin
    assert(isa(obj, 'ds.dynamicalSystem'), 'first argument must be a valid dynamicalSystems object');
    n         = obj.d.x;
    d         = obj.d.y;
    emiParams = obj.par.emiNLParams;
    assert(isfield(emiParams, 'C'), 'parameter C not found in emission parameters');
    if ~isfield(emiParams, 'bias'); emiParams.bias = []; end
    assert(isnumeric(emiParams.eta) && all(size(emiParams.eta) == [d, 10]), ...
        'emiParams.eta must be a matrix of size (d.y, 10)');
    assert(isnumeric(emiParams.bias) && (isempty(emiParams.bias) || all(size(emiParams.bias) == [d, 1])), ...
        'emiParams.bias must be a matrix of size (d.y, 1)');
    
    if nargin <= 2 || isempty(utpar) || isempty(fieldnames(utpar))
        alpha = 1;
        beta  = 0;
        kappa = 0;
    else
        assert(utils.is.numscal(utpar.alpha), 'alpha must be a numeric scalar');
        assert(utils.is.numscal(utpar.beta), 'beta must be a numeric scalar');
        assert(utils.is.numscal(utpar.kappa), 'kappa must be a numeric scalar');
        alpha = utpar.alpha; beta = utpar.beta; kappa = utpar.kappa;
    end
    
    optsDefault.etaMask    = false(0);
    optsDefault.gradient   = true;
    opts                   = utils.base.processVarargin(varargin, optsDefault);
    
    % get number of knots (excl masked-out ones)
    l            = numel(emiParams.bspl.t)-2;
    if isempty(opts.etaMask); opts.etaMask = true(1, l); end
    assert(numel(opts.etaMask)==l, 'etaMask is not same length as spline basis (%d)', l);
    keepIdx      = opts.etaMask(:);
    keepIdx2D    = repmat(keepIdx, d, 1);
    lAdapt       = sum(keepIdx);
    
    % --------- DEBUG OPTIONS ---------------------------------------
    doSlow       = false;   % for debugging
    doConstant   = false;    % constant terms in function value (ie. y*y')
    
    %% Update Parameters
    assert(all(size(x) == [d*lAdapt + d*n + d,1]), 'input must be %d by %d', [d*lAdapt + d*n + d,1]);
    
    % save current parameters
    prevEta   = obj.par.emiNLParams.eta;
    prevC     = obj.par.emiNLParams.C;
    prevBias  = obj.par.emiNLParams.bias;
    
    ds.utilIONLDS.updateParams_bspl(obj, x, keepIdx');

    emiParams = obj.par.emiNLParams;
    
    %% Calculations
    
    % SIGMA POINTS
    [W, spts]    = getSigmaPtStuff(obj, alpha, beta, kappa);
    
    % pre-calculations (+ missing vals)
    y          = obj.y;
    ymiss      = isnan(y);
    y(ymiss)   = 0;                      % non-observed (whether partial or full) y's do not contribute!
    Rinv       = inv(obj.par.R);
    RinvY      = Rinv*y;
    
    % pre-allocation
    v          = zeros(l*d, 1);
    tmpK       = zeros(l*d, (2*n + 1)*obj.d.T);   % pretty big....
    wc         = W.c;
    Hx         = zeros(d, (2*n + 1)*obj.d.T);
    
    %% MAIN LOOP (over dimension of y <- most efficient)
    %   ----------- (Major work of function) ----------------
    for dd = 1:d
        basis                    = emiParams.bspl.basisEval(spts.CXSP(dd,:));
        basis(repelem(ymiss(dd,:), 1, 2*n+1),:) = 0;   % zero missing vals (because of R^-1, cannot just do in y_t).
        
        Hx(dd,:)                 = basis * emiParams.eta(dd,:)';
        
        % basis is (2*n+1)*T   X   ell          dimension matrix
        bw                       = bsxfun(@times, basis, repmat(wc, obj.d.T, 1));
        
        % fill in relevant dimension block in v
        v((dd-1)*l+1:dd*l)       = (repelem(RinvY(dd,:), 1, 2*n+1) * bw)';
        
        % save basis elements for calculation of K (outer prod)
        tmpK((dd-1)*l+1:dd*l,:)  = basis';
    end
    K            = bsxfun(@times, tmpK, repmat(wc', 1, obj.d.T)) * tmpK';    % sum w tmpK * tmpK'
    K            = K .* kron(Rinv, ones(l));
    
    %%
    % ________________ GRADIENT FOR ETA __________________________________
    % Use quantities to calculate gradient
    thetavec     = reshape(emiParams.eta', [], 1);
    grad_eta     = v' - thetavec'*K;
    grad_eta     = grad_eta(keepIdx2D);
    
    grad_eta     = reshape(grad_eta,sum(keepIdx), d)';    % reorder
    grad_eta     = grad_eta(:);                           % reorder
    
    % ________________ FUNCTION VAL (REUSE ETA GRADIENT) __________________
    % function value
    yy           = 0;
    if doConstant || doSlow
        for tt = 1:obj.d.T
            yy = yy + y(:,tt)'*Rinv*y(:,tt);
        end
    end
    f            = v'*thetavec - 0.5*thetavec'*K*thetavec - 0.5*yy - 0.5*sum(ymiss(:));
    
    % ________________ GRADIENT FOR C, beta _______________________________
    if opts.gradient
        nspts        = 2*n+1;
        Db           = derivativeBspline(spts.CXSP, emiParams.bspl, emiParams.eta);
        ymHx         = repelem(y, 1, nspts) - Hx;
        tmp          = Db .* (Rinv*ymHx);
        wXSP         = bsxfun(@times, [spts.XSP; ones(1, nspts*obj.d.T)], repmat(wc', 1, obj.d.T));

        grad_Cb      = tmp * wXSP';
    end
    
    %% --- DEBUG -------------- (SLOW ITERATIVE EVALUATION sim. to latex)
    if doSlow
        % gradient
        [sv,sK] = slowGrad(obj, spts.CXSP, wc);
        sgrad   = sv' - thetavec'*sK;
        if norm(grad_eta - sgrad) > 1e-8
            warning('gradients differ by %e (L2 norm)', norm(grad_eta - sgrad))
        end
        
        % function val
        ff = 0;
        nspts = (2*n + 1);
        for tt = 1:obj.d.T
            curMiss = ymiss(:,tt);
            if true %~all(curMiss)
                cRinv   = Rinv;
                cRinv(curMiss,:) = [];
                cRinv(:,curMiss) = [];
                yt               = y(:,tt);
                yt(curMiss)      = [];
                for ii = 1:nspts
                    Yit    = emiParams.bspl.functionEval(spts.CXSP(:,(tt-1)*nspts+ii), emiParams.eta);
                    Yit(ymiss(:,tt)) = [];
                    ff = ff - 0.5* wc(ii) * (yt - Yit)'*cRinv*(yt - Yit);
                    ff = ff - 0.5* wc(ii) * sum(curMiss);
                end
            end
        end
        
        if norm(f - ff) > 1e-8
            warning('fun val differs by %e', norm(f - ff))
        end
        
        f    = ff;
        grad_eta = sgrad;
        
        sgrad_Cb     = 0;
        for tt = 1:obj.d.T
            for ii = 1:nspts
                ix       = (tt-1)*nspts+ii;
                delta    = y(:,tt) - Hx(:,ix);
                delta(ymiss(:,tt)) = 0;
                sgrad_Cb = sgrad_Cb + wc(ii) * (Db(:,ix) .* (Rinv*delta)) * [spts.XSP(:,ix); 1]';
            end
        end

        if norm(grad_Cb - sgrad_Cb) > 1e-8
            warning('Cb gradients differ by %e (L2 norm)', norm(grad_Cb - sgrad_Cb))
        end
        grad_Cb = sgrad_Cb;
    end
    %--- /END (DEBUG) --------------------------------------------------
    
    %% Output
    
    % restore previous values
    obj.par.emiNLParams.eta  = prevEta;
    obj.par.emiNLParams.C    = prevC;
    obj.par.emiNLParams.bias = prevBias;
    
    % return values
    if opts.gradient
        grad   = struct('eta', grad_eta, 'C', grad_Cb(:,1:n), 'bias', grad_Cb(:,n+1));
    end
    if nargout > 2
        more.v  = v(keepIdx2D);
        more.K  = K(keepIdx2D, keepIdx2D);
        more.Hx = Hx;
    end
%     fCXSP        = bspl.functionEval(spts.CXSP, eta);  
    
end


function D = derivativeBspline(X, bspl, eta)
        deriv = 1;
        ltK          = X < bspl.t(1);
        gtK          = X > bspl.t(end);
        X(ltK)       = bspl.t(1);
        X(gtK)       = bspl.t(end);
        N            = size(X,1);
        D            = zeros(size(X));
        knots        = bspl.t;
        % fda package will add in boundary knots, so remove pre-existing ones.
        knots(1:(bspl.bdknots-1)) = [];
        knots((end-(bspl.bdknots-1)+1):end) = [];
        % use fda package (Ramsay) to calculate gradient. (No time to understand it at present).
        for ii = 1:N
            dB      = external.fda.bsplineM(X(ii,:), knots, bspl.k+1, deriv, 0); 
            D(ii,:) = dB*eta(ii,:)';
        end
end

function [W, sigm] = getSigmaPtStuff(obj, alpha, beta, kappa)
    % Get sigma points sigm.XSP,
    %     weights      W.m, W.c
    %     Cx + c       sigm.CXSP
    %     Y (=h(CX+c)) sigm.Hx
    %     y - h(x)     sigm.ymHx
    
    n         = obj.d.x;
    d         = obj.d.y;
    emiParams = obj.par.emiNLParams;
    
    hasBias   = ~isempty(emiParams.bias);
    
    %% get sigma point stuff
    
    lambda   = alpha^2 * (n + kappa) - n;
    scl      = sqrt(n + lambda);
    

    % ________SIGMA POINTS X_{it}_______________________
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
        XSP(:,(1+2*n)*(tt-1) + (1:2*n+1)) = spts;
    end
    % ___________________________________________________
    
    % argument (Cx + b)
    CXSP   = obj.par.emiNLParams.C * XSP;
    if hasBias; CXSP = bsxfun(@plus, CXSP, emiParams.bias); end
    
    % attach to struct
    sigm      = struct;
    sigm.XSP  = XSP;
    sigm.CXSP = CXSP;
%     sigm.Hx   = Hx;
%     sigm.ymHx = repelem(obj.y, 1, 2*n+1) - sigm.Hx;
    
    % weights
    W         = struct;
    W.m       = ones(2*n+1, 1)./(2*(n + lambda));
    W.m(1)    = lambda/(n+lambda);
    W.c       = W.m;
    W.c(1)    = W.c(1) + (1 - alpha^2 + beta);
end

function [v,K] = slowGrad(obj, CXSP, w)
    bspl  = obj.par.emiNLParams.bspl;
    l     = numel(bspl.t)-2;
    nspts = 2*obj.d.x+1;
    d     = obj.d.y;
    
    Rinv  = inv(obj.par.R);
    ymiss = isnan(obj.y);
    y     = obj.y; y(ymiss) = 0;
    
    v     = zeros(d*l, 1);
    K     = zeros(d*l, d*l);
    for tt = 1:obj.d.T
        for ii = 1:nspts
            H = zeros(d, d*l);
            for dd = 1:d
                if ~ymiss(dd,tt)
                    B          = bspl.basisEval(CXSP(dd,(tt-1)*nspts+ii));
                    H(dd,(dd-1)*l+1:dd*l) = B;
                end
            end
            v    = v + w(ii)*(H'*Rinv*y(:,tt));
            K    = K + w(ii)*H'*Rinv*H;
        end
    end
end
