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
    end
    
    % gradient of Q
    if isempty(par) || any(strcmpi(par, 'Q'))
        A    = obj.par.A;
        D.Q  = 0.5*T * Qinv * (-1 + (s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A')*Qinv);
    end
    
    % gradient of H
    if isempty(par) || any(strcmpi(par, 'H'))
        D.H  = T * Rinv * (s.B - obj.par.H * s.SIGMA);
    end
    
    % gradient of R
    if isempty(par) || any(strcmpi(par, 'R'))
        H    = obj.par.H;
        D.R  = 0.5*T * Rinv * (-1 + (s.D - s.B*H' - H*s.B' + H*s.SIGMA*H')*Rinv);
    end
    
    % check against numerical gradient
    if doCheck
        checkGrad(obj, 'A', D.A);
        checkGrad(obj, 'Q', D.Q);
        checkGrad(obj, 'H', D.H);
        checkGrad(obj, 'R', D.R);
    end
end



function checkGrad(dso, par, exactD)
    checkGradQ(dso, par, exactD);
end

function checkGradL(dso, par, exactD)
    llhCalc  = @(x) numgradLLH(dso, x, par);
    approxD  = utils.math.numgrad(llhCalc, dso.par.(par), 1e-5);
    fprintf('**Param ''%s''**: ''exact'' gradient:\n', par); disp(exactD);
    fprintf('numerical gradient:\n'); disp(approxD);
    fprintf('distance (L2) = %3.6f\n', norm(exactD - approxD));
end

function checkGradQ(dso, par, exactD)
    llhCalc  = @(x) numgradQ(dso, x, par);
    approxD  = utils.math.numgrad(llhCalc, dso.par.(par), 1e-5);
    fprintf('**Param ''%s''**: ''exact'' gradient:\n', par); disp(exactD);
    fprintf('numerical gradient:\n'); disp(approxD);
    fprintf('distance (L2) = %3.6f\n', norm(exactD - approxD));
end


function llh = numgradLLH(dso, X, par)
    dso.par.(par) = X;
    dso = dso.filterKalman(true, false);
    llh = dso.infer.llh;
end

function q = numgradQ(dso, X, par)
    dso.par.(par) = X;
    dso = dso.filterKalman(false, false);
%     dso = dso.smoothLinear(false);
    q = dso.expLogJoint;
end