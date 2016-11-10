function D = getGradient(obj, par, doCheck)
    
    if nargin < 3 || isempty(doCheck)
        doCheck = false;
    end
    
    if nargin < 2 || isempty(par)
        par = {};
    end
    s    = obj.suffStats;
    T    = obj.d.T;
    
    D    = struct;
    
    % pre-compute matrix inverses
    if isempty(par) || any(strcmpi(par, 'A')) || any(strcmpi(par, 'Q'))
        Qinv = inv(obj.par.Q);
    end
    if isempty(par) || any(strcmpi(par, 'H')) || any(strcmpi(par, 'R'))
        Rinv = inv(obj.par.R);
    end
    
    % gradient of A
    if isempty(par) || any(strcmpi(par, 'A'))
        D.A  = T * Qinv * (s.C - obj.par.A * s.PHI);
        if doCheck > 0
            checkGrad(obj, 'A', D.A, doCheck-1);
        end
    end
    
    % gradient of Q
    if isempty(par) || any(strcmpi(par, 'Q'))
        A    = obj.par.A;
        D.Q  = -0.5*T*Qinv + 0.5*T*Qinv*(s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A')*Qinv;
        if doCheck > 0
            checkGrad(obj, 'Q', D.Q, doCheck-1);
        end
    end
    
    % gradient of H
    if isempty(par) || any(strcmpi(par, 'H'))
        D.H  = T * Rinv * (s.B - obj.par.H * s.SIGMA);
        if doCheck > 0
            checkGrad(obj, 'H', D.H, doCheck-1);
        end
    end
    
    % gradient of R
    if isempty(par) || any(strcmpi(par, 'R'))
        H    = obj.par.H;
        D.R  = -0.5*T*Rinv + 0.5*T*Rinv*(s.D - s.B*H' - H*s.B' + H*s.SIGMA*H')*Rinv;
        if doCheck > 0
            checkGrad(obj, 'R', D.R, doCheck-1);
        end
    end
    
%     % check against numerical gradient
%     if doCheck
%         checkGrad(obj, 'A', D.A);
%         checkGrad(obj, 'Q', D.Q);
%         checkGrad(obj, 'H', D.H);
%         checkGrad(obj, 'R', D.R);
%     end
end


function checkGrad(dso, par, exactD, verbose)
    llhCalc  = @(x) numgradQ(dso, x, par);
%     llhCalc  = @(x) numgradLLH(dso, x, par); % doesn't work so well -> must do inference on nonPSD matrix
    approxD  = utils.math.numgrad(llhCalc, dso.par.(par), 1e-5);
    nL2      = norm(exactD - approxD)./numel(exactD);
    if verbose
        fprintf('**Param ''%s''**: ''exact'' gradient:\n', par); disp(exactD);
        fprintf('numerical gradient:\n'); disp(approxD);
        fprintf('distance (L2) = %3.6f\n', nL2);
    else
        fprintf('(%s) | distance (L2) = %3.6f\n', par, nL2)
    end
end


function llh = numgradLLH(dso, X, par)
    dso.par.(par) = X;
    fOpts = struct('bDoValidation', false, 'bIgnoreHash', true);
    dso = dso.filter('Kalman', true, [], fOpts);
    dso = dso.smooth('Linear', [], fOpts);
    llh = dso.infer.llh;
end

function q = numgradQ(dso, X, par)
    dso.par.(par) = X;
    % -- Note that we must differentiate Q *evaluated at* current --
    % -- value of theta. Thus while we must perturb the parameters --
    % -- to evaluate the gradient in the log joint, we don't change --
    % -- the posterior. We must hack the posterior so it doesn't --
    % -- notice the change! --
    dso.infer.fpHash = dso.parameterHash;  % <-- Hack!
    q = dso.expLogJoint;
end