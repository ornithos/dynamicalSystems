function [obj, llh] = parameterLearningEM(obj, opts)

optsDefault     = struct('epsilon', 1e-3, 'maxiter', 200, 'ssid', 0, 'verbose', true);
opts            = utils.base.parse_argumentlist(opts, optsDefault);

% Set initial values
if opts.ssid
    fprintf('(%s) Initialising using subspace identification...\n', datestr(now, 'HH:MM:SS'));
    obj = ssid(obj, opts.ssid);
else
    noWarningNeeded = false;
    var_y   = var(obj.y);
    if isempty(obj.A)
        if ~noWarningNeeded && ~isempty(obj.inA)
            warning('initialising with null A rather than value given in construction.'); 
            noWarningNeeded = true; 
        end
        obj.A = eye(obj.d.x);
    end
    if isempty(obj.Q)
        if ~noWarningNeeded && ~isempty(obj.inQ)
            warning('initialising with null Q rather than value given in construction.'); 
            noWarningNeeded = true; 
        end
        obj.Q = var_y*eye(obj.d.x)/10;
    end
    if isempty(obj.H)
        if ~noWarningNeeded && ~isempty(obj.inH)
            warning('initialising with null H rather than value given in construction.'); 
            noWarningNeeded = true; 
        end
        obj.H = eye(obj.d.x);
    end
    if isempty(obj.R)
        if ~noWarningNeeded && ~isempty(obj.inR)
            warning('initialising with null R rather than value given in construction.');
        end
        obj.R = eye(obj.d.y);
    end
end

llh        = [-Inf; zeros(opts.maxiter,1)];
% prevParams = getAllParams(obj);

for ii = 1:opts.maxiter
    obj    = obj.filterKalman(true, false);
    obj    = obj.smoothLinear;

    llh(ii+1) = obj.llh;
    delta     = llh(ii+1) - llh(ii);
    if abs(delta) < opts.epsilon
        fprintf('(%s) EM Converged in %d iterations (%.4f) \n', datestr(now), ii, delta);
        break
    end
%     obj       = obj.parameterLearningMStep;
    obj       = obj.parameterLearningMStep([], {'A', 'Q'});
    obj       = obj.filterKalman;
    obj       = obj.smoothLinear;
    obj       = obj.parameterLearningMStep([], {'H', 'R'});
    
    if opts.verbose
        fprintf('--- (%s) Iteration %4d: LLH change: % .2f\n', ...
            datestr(now, 'HH:MM:SS'), ii, delta);
    end
end
llh = llh(2:ii+1);
end

function out = getAllParams(s)
    out = struct;
    out.A = s.A;
    out.H = s.H;
    out.Q = s.Q;
    out.R = s.R;
end

function [out, which] = testParams(q1, q2)
    mse = zeros(4,1);
    strTheta = {'A','H','Q','R'};
    mse(1) = sum(sum((q1.A - q2.A).^2));
    mse(2) = sum(sum((q1.H - q2.H).^2));
    mse(3) = sum(sum((q1.Q - q2.Q).^2));
    mse(4) = sum(sum((q1.R - q2.R).^2));
    out    = sum(sqrt(mse));
    [~,which] = max(mse);
    which  = strTheta{which};
end