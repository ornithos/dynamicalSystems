function [obj, llh, ii] = parameterLearningEM(obj, opts)

if nargin < 2 || isempty(opts); opts = struct; end
optsDefault     = struct('epsilon', 1e-3, 'maxiter', 200, 'ssid', false, 'ssidL', 5, ...
                        'verbose', true, 'dbg', false, 'validation', false, ...
                        'sampleStability', 1);
optsDefault     = utils.base.parse_argumentlist(obj.opts, optsDefault, false);      % bring in global opts
opts            = utils.base.parse_argumentlist(opts, optsDefault, false);          % add user specified opts.

sing_val_eps    = 0.004999999999;  % tolerance of singular values of A > 1. This value
                                   % is given in in Sidiqqi et als code for constraint
                                   % generation.

% Set initial values
if opts.ssid
    if opts.verbose
        fprintf('(%s) Initialising using subspace identification...\n', datestr(now, 'HH:MM:SS'));
    end
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

if opts.validation
    obj.validationInference;  % ensure able to do inference
end

fOpts      = struct('bDoValidation', false, 'bIgnoreHash', true);
converged  = false;
iterBar    = utils.base.objIterationBar;
iterBar.newIterationBar('EM Iteration: ', opts.maxiter, true, '--- ', 'LLH change: ');


% multistep: do all 4 maximisations together initially. Initialise dbg struct in case we need it..
multiStep  = 4;
dbgLLH = struct('A',[0,0],'Q',[0,0],'H',[0,0],'R',[0,-Inf]);
for ii = 1:opts.maxiter
    obj    = obj.filter('Kalman', true, [], fOpts);
    obj    = obj.smooth('Linear', [], fOpts);
    
    % llh calc
    dbgLLH.R  = [obj.infer.llh - dbgLLH.H(2), obj.infer.llh];
    llh(ii+1) = obj.infer.llh;
    delta     = llh(ii+1) - llh(ii);
    
    % negative update warning
    if delta < -1e-8
        iterBar.updateText([iterBar.text, '*']);
        if multiStep == 1
            % basically things have gone really wrong by here..
            iterBar.clearConsole;
            fprintf('min eigv Q = %.3e', min(eig(obj.par.Q)));
            fprintf(', min eigv R = %.3e\n', min(eig(obj.par.R)));
            keyboard
        elseif opts.sampleStability > 1
            % only periodically sampling stability of A leads to jumps..
            opts.sampleStability = 1;
        else
            % Doing ECM is not *guaranteed* to increase the llh on each step
            multiStep = multiStep ./ 2;
        end
    end
    
    % convergence
    if abs(delta) < opts.epsilon
        iterBar.finish;
        converged = true;
        if opts.verbose
            fprintf('(%s) EM Converged in %d iterations (%.3e < %.3e) \n', datestr(now), ii, delta, opts.epsilon);
        end
        break
    end
    
    % M-Steps (interleaved with E-steps if required, since ECM not EM!)
    % ----------------------------
    
    % ____ Canonical parameters ___________________________________________
    obj       = obj.parameterLearningMStep(opts.dbg, {'A'});
    
    % Check for stability of A: significant problems when A blows up.
    if mod(ii, opts.sampleStability) == 0
        if max(abs(eig(obj.par.A))) > 1 + sing_val_eps
            cgTic        = tic;
            if opts.verbose
                iterBar.clearConsole;
                fprintf('Stabilising A with constraint generation... ');
            end
            obj.par.A    = stabiliseA_constraintGeneration(obj, opts.verbose); % see mini function below ????
            if opts.verbose; fprintf('Done! (%.2f)\n', toc(cgTic)); end
        end
    end
    
    % Q, H, R
    switch multiStep
        case 4
            obj       = obj.parameterLearningMStep(opts.dbg, {'Q','H','R'});
        case 2
            obj       = obj.parameterLearningMStep(opts.dbg, {'Q'});
            obj       = obj.filter('Kalman', true, [], fOpts);
            obj       = obj.smooth('Linear', [], fOpts);
            obj       = obj.parameterLearningMStep(opts.dbg, {'H','R'});
        case 1
            obj       = obj.filter('Kalman', true, [], fOpts);
            obj       = obj.smooth('Linear', [], fOpts);
            dbgLLH.A  = [obj.infer.llh - dbgLLH.R(2), obj.infer.llh];
            obj       = obj.parameterLearningMStep(opts.dbg, {'Q'});
            obj       = obj.filter('Kalman', true, [], fOpts);
            obj       = obj.smooth('Linear', [], fOpts);
            dbgLLH.Q  = [obj.infer.llh - dbgLLH.A(2), obj.infer.llh];
            obj       = obj.parameterLearningMStep(opts.dbg, {'H'});
            obj       = obj.filter('Kalman', true, [], fOpts);
            obj       = obj.smooth('Linear', [], fOpts);
            dbgLLH.H  = [obj.infer.llh - dbgLLH.Q(2), obj.infer.llh];
            obj       = obj.parameterLearningMStep(opts.dbg, {'R'});
        otherwise
            error('Param Learning: FAIL. multiStep NOT IN (1,2,4)');
    end
    
    % ____ Control parameters _____________________________________________
    if any(obj.hasControl)
        obj       = obj.filter('Kalman', [], [], fOpts);
        obj       = obj.smooth('Linear', [], fOpts);
        dbgLLH.R  = [obj.infer.llh - dbgLLH.H(2), obj.infer.llh];
        if multiStep > 1
            obj       = obj.parameterLearningMStep(opts.dbg, {'B', 'C'});
        else
            obj       = obj.parameterLearningMStep(opts.dbg, {'B'});
            obj       = obj.filter('Kalman', [], [], fOpts);
            obj       = obj.smooth('Linear', [], fOpts);
            dbgLLH.B  = [obj.infer.llh - dbgLLH.Q(2), obj.infer.llh];
            dbgLLH.H  = [dbgLLH.B(1), dbgLLH.H(2)];
            obj       = obj.parameterLearningMStep(opts.dbg, {'C'});
        end
    end
    
    % update console
    if opts.verbose && ~opts.dbg
        iterBar.print(ii, delta);
    end
end

if ~converged && opts.verbose
    iterBar.finish;
    fprintf('EM did not converge in %d iterations\n', opts.maxiter);
end

llh = llh(2:ii+1);
end


function A = stabiliseA_constraintGeneration(obj, verbose)
    % -- what about x0 estimates?
    m          = obj.infer.smooth.mu;
    P          = obj.infer.smooth.sigma;
    G          = obj.infer.smooth.G;
    sumP       = sum(cat(3,P{1:obj.d.T}),3);
    sumC       = 0;
    for tt = 1:obj.d.T-1
        sumC = sumC + P{tt+1} * G{tt}';
    end    

    S1        = m(:,1:end-1);
    S2        = m(:,2:end);
    sumP      = (sumP + sumP)./2;
    A         = ds.utils.learnCGModel_EM(S1, S2, sumP, sumC, obj.par.A, 0, verbose);
end
