function [obj, llh] = parameterLearningEM(obj, opts)

if nargin < 2 || isempty(opts); opts = struct; end
optsDefault     = struct('epsilon', 1e-3, 'maxiter', 200, 'ssid', false, 'ssidL', 5, 'verbose', true);
opts            = utils.base.parse_argumentlist(opts, optsDefault);

% Set initial values
if opts.ssid
    fprintf('(%s) Initialising using subspace identification...\n', datestr(now, 'HH:MM:SS'));
    obj = obj.ssid(opts.ssidL);
else
    var_y   = var(obj.y);
    if isempty(obj.par.A)
        obj.par.A = eye(obj.d.x);
    end
    if isempty(obj.par.Q)
        obj.par.Q = var_y*eye(obj.d.x)/10;
    end
    if isempty(obj.par.H)
        obj.par.H = eye(obj.d.x);
    end
    if isempty(obj.par.R)
        obj.par.R = eye(obj.d.y);
    end
end

llh        = [-Inf; zeros(opts.maxiter,1)];

obj.validationInference;  % ensure able to do inference
fOpts      = struct('bDoValidation', false, 'bIgnoreHash', true);
for ii = 1:opts.maxiter
    obj    = obj.filter('Kalman', true, [], fOpts);
    obj    = obj.smooth('Linear', [], fOpts);

    llh(ii+1) = obj.infer.llh;
    delta     = llh(ii+1) - llh(ii);
    if abs(delta) < opts.epsilon
        fprintf('(%s) EM Converged in %d iterations (%.3e < %.3e) \n', datestr(now), ii, delta, opts.epsilon);
        break
    end
%     obj       = obj.parameterLearningMStep;
    obj       = obj.parameterLearningMStep([], {'A', 'Q'});
    obj       = obj.filter('Kalman', [], [], fOpts);
    obj       = obj.smooth('Linear', [], fOpts);
    obj       = obj.parameterLearningMStep([], {'H', 'R'});
    
    if opts.verbose
        fprintf('--- (%s) Iteration %4d: LLH change: % .2f\n', ...
            datestr(now, 'HH:MM:SS'), ii, delta);
    end
end
llh = llh(2:ii+1);
end