function [a, D, q] = expLogJoint_bspl(obj, varargin)

    % expLogJoint_bspl
    % Method of an IONLDS systems object. Calculates the expected log
    % joint **given the current posterior**, <log(p(x,z))>
    % OF SPLINE NONLINEARITY!
    %
    % Usage:
    %  val           = obj.ExpLogJoint
    %          - Outputs the value of <log(p(x,z))>.
    %
    %  [val, D]      = obj.ExpLogJoint
    %          - second output gives partial derivatives wrt all params.
    %
    %  [val, D, decomp] = obj.ExpLogJoint
    %          - third output gives decomposition of value as normalising
    %            constants (1), (2), and quadratic forms of both distns as
    %            (3), (4).
    %
    %  val           = obj.expLogJoint(...,'utpar', s)
    %          - struct specifying UT parameters *alpha *beta *kappa
    %
    %  val           = obj.ExpLogJoint(...,'freeEnergy', true, ...)
    %          - Outputs the free energy <log(p(x,z))> + H(q).
    %
    %  val           = obj.ExpLogJoint(..., 'recalcPosterior', true, ...)
    %          - Recalculates the posterior before calculating expLogJoint.
    %            Note that the new posterior will not be saved to the ds
    %            object, ie. it will not affect the current state of
    %            inference.
    %
    
    optsDefault = struct('freeEnergy', false, 'recalcPosterior', false, ...
                         'utpar', [], 'gamma', false);
                         
    opts        = utils.base.processVarargin(varargin, optsDefault);
    ssopts      = struct('bIgnoreHash', ~opts.recalcPosterior);
    
    %%
    tmpobj  = obj.copy;
    s       = tmpobj.suffStats(ssopts);
    s2      = s.emissions;
    A       = tmpobj.par.A;
    %H       = tmpobj.par.H;

    q       = NaN(4,1);
    q(1)    = -0.5*s.T   * utils.math.logdet(2*pi*tmpobj.par.Q);
    q(2)    = -0.5*s2.Ty * utils.math.logdet(2*pi*tmpobj.par.R);
    
%     hasBias = ~isempty(obj.par.c);   % bias not relevant for ionlds.
    
    if ~obj.hasControl(1)
        ctrlAdd = zeros(obj.d.x);
    else
        Bum     = obj.par.B * s.XU';
        BuAm_m  = obj.par.B * s.Xm1_U' * obj.par.A';
        ctrlAdd = obj.par.B*s.UU*obj.par.B' - Bum - Bum' + BuAm_m + BuAm_m';
    end
    
    q(3)    = -0.5*s.T*trace((tmpobj.par.Q)\(s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A' + ctrlAdd));
    
    % Nonlinear bit
    % -- IONLDS = no control in emission
    if ~isa(obj, 'ds.dynamicalSystemBatch')    
        if isempty(opts.utpar), [ymHx, M2, XSP, Wc] = ds.utilIONLDS.utTransform_ymHx_bspl(obj);
        else,                   [ymHx, M2, XSP, Wc] = ds.utilIONLDS.utTransform_ymHx_bspl(obj, utpar.alpha, utpar.beta, utpar.kappa);
        end
    else
        if isempty(opts.utpar), [ymHx, M2, XSP, Wc] = ds.utilIONLDS.utTransform_ymHx_bspl_batch(obj);
        else,                   [ymHx, M2, XSP, Wc] = ds.utilIONLDS.utTransform_ymHx_bspl_batch(obj, utpar.alpha, utpar.beta, utpar.kappa);
        end
    end
        
    q(4)    = -0.5*trace((tmpobj.par.R)\M2);
    
    % x0
    parV0   = tmpobj.par.x0.sigma;
    parm0   = tmpobj.par.x0.mu;
    P0      = tmpobj.infer.smooth.x0.sigma;
    m0      = tmpobj.infer.smooth.x0.mu;
    q(5)    = -0.5*utils.math.logdet(2*pi*parV0) -0.5*trace((parV0)\(P0 + (m0-parm0)*(m0-parm0)'));
    
    if nargout > 1 && false
        % *****************************************************************
        % NOT IMPLEMENTED YET
        % *****************************************************************
        % *****************************************************************
        
        D     = struct;
        dM2   = struct;

        eta   = obj.par.emiNLParams.eta;
        m     = eta(:,1); M = eta(:,2); nu = eta(:,4); gamma = eta(:,4);
        if ~opts.gamma  % roll gamma into C -> overparameterised.
            gamma = ones(size(gamma));
        end
        
        C     = [obj.par.emiNLParams.C];
        
        % ____ get derivatives of M2 first _______________________________
        % CXSP are C * sigma points of smoothed distn. The tensor has been
        % unrolled such that CXSP(:,1:2n+1) are the sigma points of (tt=1),
        % CXSP(:,2n+2:(2n+1)*2) are the sigma points of (tt=2) etc.
        %
        %           (useful quantities #1)
        
        %XSP          = [ones(1, size(XSP,2)); XSP];  % prepend one for bias
        CXSP         = C * XSP;
        gCX          = bsxfun(@times, gamma, CXSP);
        expmgCX      = exp(-gCX);
        denom        = 1 + expmgCX;
        denomPowNu   = bsxfun(@power, denom, 1./nu);
        
        % *** FOR QUANTITIES (M, m, nu, gamma), the Jacobian is diagonal,
        % *** and we can reduce storage by concatenating diags in matrix.
        dM2.M        = 1./denomPowNu;
        dM2.m        = 1 - dM2.M;
        
        log1PEGCX    = utils.math.log1pexp(-gCX);
        dM2.nu       = bsxfun(@times, (M-m)./(nu.^2), log1PEGCX./denomPowNu);
        
        % FOR DEBUG>>>>
%         testgrad     = utils.math.numgrad(@(x) testDeriv(obj, CXSP, x, 1, 3), 0, 1e-8);
        
        %           (useful quantities #2)
        denom2       = (1+exp(gCX)).*denom;
        tmpGC        = bsxfun(@times, (M-m)./(nu), 1./denom2);
        
        dM2.gamma    = tmpGC .* CXSP;
        dM2.C        = bsxfun(@times, tmpGC, gamma);   

        % ______ Get sigma point weighted sum ____________________________
        Wc           = repmat(Wc', 1, T);   % matches unrolled CXSP
        WymHx        = bsxfun(@times, Wc, ymHx); 
        
        % FOR DEBUG>>>>
%         testgrad     = utils.math.numgrad(@(x) testDerivM(obj, CXSP, Wc(1:5)', x, 1, 4), 0, 1e-8);
                
        dM2.m        = sum(WymHx .* dM2.m, 2);
        dM2.M        = sum(WymHx .* dM2.M, 2);
        dM2.nu       = sum(WymHx .* dM2.nu, 2);
        dM2.gamma    = sum(WymHx .* dM2.gamma, 2);

        % ______ derivative = -0.5 * (-2) * R^{-1} * d/dpsi M2 ___________
        Dtmp         = (tmpobj.par.R)\ [dM2.m, dM2.M, dM2.nu, dM2.gamma];
        D.m          = Dtmp(:,1);
        D.M          = Dtmp(:,2);
        D.nu         = Dtmp(:,3);
        D.gamma      = Dtmp(:,4);
        %if ~opts.gamma; D.gamma = ones(size(D.gamma)); end

        % ================================================================
        % ================================================================
        % *** For C, the Jacobian is 3 dimensional,
        % *** so we have to take a different approach.
        
        D.C         = (WymHx .* dM2.C) * XSP';
        D.C         = (tmpobj.par.R) \ D.C;
        
        D.b         = D.C(1,:);
        D.C         = D.C(2:end,:);
        % probably could have done this for them all tbh 
    else
        D = [];
    end
    
    if opts.freeEnergy
        % entropy calculation for free energy
        if obj.d.T*obj.d.x > 2000
            warning('Determinant of size %d x %d will be calculated...', obj.d.T*obj.d.x,obj.d.T*obj.d.x);
        end
        fullcov = ds.utils.fullJointCovariance(obj, s.T);
        detArg  = 2*pi*fullcov;
        L       = chol(detArg);
        q(6)    = 0.5*2*sum(log(diag(L))) + (s.T + 1)*obj.d.x./2;  % logdet = 2*sum(log(diag(chol(.))))
        % T + 1 since posterior from 0:T.
        
        % -------------------------------------------------
        % (MISSING VALUES!!)entropy of distribution over y's 
        nans    = isnan(obj.y);
        if true && sum(sum(nans)) > 0
            q(7)   = 0;
            R      = obj.par.R;
            nanIdx = find(sum(nans)>0 & ~all(nans));
            for ii = 1:numel(nanIdx)
                jj     = nanIdx(ii);
                Rtmp   = R(nans(:,jj), nans(:,jj));
                q(7)   = q(7) + 0.5*utils.math.logdet(2*pi*Rtmp*exp(1));
            end
        end
    end
    
    a       = sum(q);
end


function out = testDeriv(obj, CXSP, test, iy, ix)
    eta     = obj.par.emiNLParams.eta;
    
    eta(iy, ix) = eta(iy, ix) + test;
    out     = ds.utilIONLDS.gen_sigmoid(CXSP, eta);
end


function out = testDerivM(obj, CXSP, Wc, test, iy, ix)
    eta     = obj.par.emiNLParams.eta;
    eta(iy, ix) = eta(iy, ix) + test;
    y       = repelem(obj.y, 1, numel(Wc));
    Wc      = repmat(Wc', 1, obj.d.T);
    out     = sum(Wc .* sum((y - ds.utilIONLDS.gen_sigmoid(CXSP, eta)).^2, 1),2);
end

function out = allEtaGrad(fn, ny, nx)
    out   = zeros(ny, nx);
    for ii = 1:nx
        for jj = 1:ny
            out(jj,ii) = fn(jj, ii);
        end
    end
end