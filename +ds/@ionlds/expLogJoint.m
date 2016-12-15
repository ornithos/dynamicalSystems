function [a, D] = expLogJoint(obj, utpar)

    if nargin < 2
        utpar = [];
    end

    tmpobj = obj;
    opts = struct('bIgnoreHash', true);
    s    = tmpobj.suffStats(opts);
    A    = tmpobj.par.A;
    H    = tmpobj.par.H;
    T    = tmpobj.d.T;
    q    = NaN(4,1);
    q(1)    = -0.5*T*log(det(2*pi*tmpobj.par.Q));
    q(2)    = -0.5*T*log(det(2*pi*tmpobj.par.R));
    
    if ~obj.hasControl(1)
        ctrlAdd = zeros(obj.d.x);
    else
        Bum     = obj.par.B * s.XU';
        BuAm_m  = obj.par.B * s.Xm1_U' * obj.par.A';
        ctrlAdd = obj.par.B*s.UU*obj.par.B' - Bum - Bum' - BuAm_m - BuAm_m';
    end
    
    q(3)    = -0.5*T*trace((tmpobj.par.Q)\(s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A' + ctrlAdd));
    
    % Nonlinear bit
    % -- IONLDS = no control in emission
    if isempty(utpar), [ymHx, M2, XSP, Wc] = ds.utilIONLDS.utTransform_ymHx(obj);
    else,              [ymHx, M2, XSP, Wc] = ds.utilIONLDS.utTransform_ymHx(obj, utpar.alpha, utpar.beta, utpar.kappa);
    end
        
    q(4)    = -0.5*trace((tmpobj.par.R)\M2);
    a       = sum(q);
    
    if nargout > 1
        D     = struct;
        dM2   = struct;

        eta   = obj.par.emiNLParams.eta;
        m     = eta(:,1); M = eta(:,2); nu = eta(:,3); gamma = eta(:,4);
        b     = eta(:,5);
        C     = [b, obj.par.emiNLParams.C];
        
        % ____ get derivatives of M2 first _______________________________
        % CXSP are C * sigma points of smoothed distn. The tensor has been
        % unrolled such that CXSP(:,1:2n+1) are the sigma points of (tt=1),
        % CXSP(:,2n+2:(2n+1)*2) are the sigma points of (tt=2) etc.
        %
        %           (useful quantities #1)
        
        XSP          = [ones(1, size(XSP,2)); XSP];  % prepend one for bias
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

        % ================================================================
        % ================================================================
        % *** For C, the Jacobian is 3 dimensional,
        % *** so we have to take a different approach.
        
        D.C         = (WymHx .* dM2.C) * XSP';
        D.C         = (tmpobj.par.R) \ D.C;
        
        D.b         = D.C(1,:);
        D.C         = D.C(2:end,:);
        % probably could have done this for them all tbh 
    end
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