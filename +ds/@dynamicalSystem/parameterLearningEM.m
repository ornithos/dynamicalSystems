function [llh, ii] = parameterLearningEM(obj, opts)

if nargin < 2 || isempty(opts); opts = struct; end
optsDefault     = struct('epsilon', 1e-3, 'maxiter', 200, 'ssid', false, 'ssidL', 5, ...
                        'verbose', true, 'dbg', false, 'validation', false, ...
                        'multistep', 4, 'diagQ', false, 'diagR', false, ...
                        'sampleStability', 1, 'stableVerbose', false, ...
                        'annealingSchedule', Inf, 'annealingIter', 10, ...
                        'annealingMin', 1e-6, 'strictNegativeCheck', false);
optsDefault     = utils.base.parse_argumentlist(obj.opts, optsDefault, false);      % bring in global opts
opts            = utils.base.parse_argumentlist(opts, optsDefault, false);          % add user specified opts.

sing_val_eps    = 0.004999999999;  % tolerance of singular values of A > 1. This value
                                   % is given in in Sidiqqi et als code for constraint
                                   % generation.

% ________ Set initial values ____________________________________________
if opts.ssid
    if opts.verbose
        fprintf('(%s) Initialising using subspace identification...\n', datestr(now, 'HH:MM:SS'));
    end
    obj.ssid(opts.ssidL);
    ssidrescale = max(abs(eig(obj.par.A)));
    ssidrescale = max(ssidrescale, 0);    % ensure stable trans matrix (hack!)
    obj.par.A   = obj.par.A ./ ssidrescale;
    obj.par.H   = obj.par.H .* ssidrescale;
    
    obj.filter;
    obj.smooth;
%     % do filter and smooth first anyway, since need suff stats for constrain gen
%     if max(abs(eig(obj.par.A))) > 1 + sing_val_eps
%         obj.par.A    = stabiliseA_constraintGeneration(obj, obj.par.A, 1);
%         obj.filter;
%         obj.smooth;
%     end
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

if opts.validation
    obj.validationInference;  % ensure able to do inference
end

% ________ Determine annealing schedule __________________________________
% Perform Deterministic Annealing (cf. Ueda, Nakano - Deterministic
% Annealing Variant of the EM Algorithm, NIPS 1995)
% beta denotes the inverse temperature: 0 = hot, 1 = cooled.
assert(opts.annealingSchedule > 0, 'annealing schedule must increase by strictly positive amount');
assert(opts.annealingMin > 0 && opts.annealingMin <= 1, 'annealing minimum beta must be in (0,1]');
da.invbetas     = (linspace(1, 0, max(1,1./opts.annealingSchedule)));
da.invbetas     = exp(-da.invbetas.*log(opts.annealingMin));
da.invbetas     = round(da.invbetas,5);
da.n            = numel(da.invbetas);
da.pointer      = 1;
da.cur          = da.invbetas(1);  % may not be min if annealingSchedule is infinite.
da.curIter      = 0;

% ________ Get names of parameters which have been 'fixed' _______________
% get options for Filtering and for M-step
fOpts      = struct('bDoValidation', false, 'bIgnoreHash', true);
mstepOpts  = struct('verbose', opts.dbg, 'diagQ', opts.diagQ, 'diagR', opts.diagR);
optFds     = fieldnames(opts);
for oname = optFds   % fixA, fixQ, fix....
    if strcmp(oname(1:3), 'fix')
        mstepOpts.(oname) = opts.(oname);
    end
end

% remove linear updates where non-linear function
if ~obj.evoLinear
    mstepOpts.fixA = true;
    mstepOpts.fixB = true;
    mstepOpts.fixQ = true;
end
if ~obj.emiLinear
    mstepOpts.fixH = true;
    mstepOpts.fixC = true;
    mstepOpts.fixR = true;
end

% preallocate deterministic annealing amount
mstepOpts.anneal = 1;   % 1 corresponds to no annealing.

% _______ initialise _____________________________________________________
llh        = [-Inf; zeros(opts.maxiter,1)];
converged  = false;
iterBar    = utils.base.objIterationBar;
iterBar.newIterationBar('EM Iteration: ', opts.maxiter, true, '--- ', 'LLH change: ');


% multistep: do all 4 (default) maximisations together initially.
multiStep  = opts.multistep;
if ~(~opts.verbose && ~opts.stableVerbose)
    % if no verbose, switch off iteration count (.) in constraint gen script.
    % else stableVerbose = 1 is the dots, stableVerbose = 2 is the full thing
    opts.stableVerbose = opts.stableVerbose + 1; 
end

% Initialise dbg struct in case we need it..
dbgLLH = struct('A',[0,0],'Q',[0,0],'H',[0,0],'R',[0,-Inf]);

%% MAIN EM LOOP
for ii = 1:opts.maxiter
    % E-Step!
    obj.filter('Kalman', true, [], fOpts);
    obj.smooth('Linear', [], fOpts);
    
    % llh calc
    dbgLLH.R  = [obj.infer.llh - dbgLLH.H(2), obj.infer.llh];
    llh(ii+1) = obj.infer.llh;
    delta     = llh(ii+1) - llh(ii);
    
    % ============== Convergence and Admin ===============================
    % negative update warning
    % --- if prev step constrained, we may have seen a drop in LLH.
    %     However, we do not expect this in the case where sampleStability
    %     is 1 since should still be monotonic if every iter stable.
    % --- Deterministic annealing artificially reflates the variance so no
    %     guarantees if beta < 1.
    if opts.strictNegativeCheck
        negativeLlhStep = dbgLLH.R(1) < -1e-8 || dbgLLH.Q(1) < -1e-8 || dbgLLH.A(1) < -1e-8 || dbgLLH.H(1) < 1e-8;
    else
        negativeLlhStep = delta < -1e-8 && da.cur == 1 && ~(prevStepWasConstrained && ~(opts.sampleStability == 1));
    end
    
    if negativeLlhStep && ii > 1
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
        % annealing has finished?
        if da.cur == 1
            if opts.verbose; iterBar.finish; end
            converged = true;
            if opts.verbose
                fprintf('(%s) EM Converged in %d iterations (%.3e < %.3e) \n', datestr(now), ii, delta, opts.epsilon);
            end
            break
        else
            % update annealing if not
            da.pointer = min(da.pointer + 1, da.n);
            da.cur     = da.invbetas(da.pointer);
            da.curIter = 1;
            iterBar.updateText([iterBar.text, '>']);
        end
    end
    % ====================================================================
    %
    %% M-Steps (interleaved with E-steps if required, since ow. ECM not EM!)
    % ----------------------------
    mstepOpts.anneal = da.cur;
    
    % ____ Canonical parameters ___________________________________________
    prevA         = obj.par.A;
%     prevLLH       = obj.logLikelihood;
    obj.parameterLearningMStep({'A', 'B'}, mstepOpts);
    
    % Check for stability of A: significant problems when A blows up.
    prevStepWasConstrained = false;
    if mod(ii, opts.sampleStability) == 0
        if max(abs(eig(obj.par.A))) > 1 + sing_val_eps
            cgTic        = tic;
            if opts.stableVerbose > 1
                iterBar.clearConsole;
                fprintf('Stabilising A with constraint generation... ');
            end
            badA         = obj.par.A;
            obj.par.A    = prevA;   % reset to last good A
            obj.par.A    = stabiliseA_constraintGeneration(obj, prevA, opts.stableVerbose); % see mini function below ????
            if max(abs(eig(obj.par.A))) > 1
                 iterBar.updateText([iterBar.text, '*']);
%                 keyboard
            end
%             newLLH       = obj.logLikelihood;
%             if newLLH < prevLLH
%                 keyboard
%             end
            if opts.stableVerbose; fprintf('Done! (%.2f)\n', toc(cgTic)); end
            prevStepWasConstrained = true;
        end
    end
    
    % Q, H, R
    switch multiStep
        case 4
            obj.parameterLearningMStep({'Q','H','R'}, mstepOpts);
        case 2
            obj.parameterLearningMStep({'Q'}, mstepOpts);
            obj.filter('Kalman', true, [], fOpts);
            obj.smooth('Linear', [], fOpts);
            obj.parameterLearningMStep({'H','R'}, mstepOpts);
        case 1
            obj.filter('Kalman', true, [], fOpts);
            obj.smooth('Linear', [], fOpts);
            dbgLLH.A  = [obj.infer.llh - dbgLLH.R(2), obj.infer.llh];
            
            obj.parameterLearningMStep({'Q'}, mstepOpts);
            obj.filter('Kalman', true, [], fOpts);
            obj.smooth('Linear', [], fOpts);
            dbgLLH.Q  = [obj.infer.llh - dbgLLH.A(2), obj.infer.llh];
            
            obj.parameterLearningMStep({'H'}, mstepOpts);
            obj.filter('Kalman', true, [], fOpts);
            obj.smooth('Linear', [], fOpts);
            dbgLLH.H  = [obj.infer.llh - dbgLLH.Q(2), obj.infer.llh];
            
            obj.parameterLearningMStep({'R'}, mstepOpts);
        otherwise
            error('Param Learning: FAIL. multiStep NOT IN (1,2,4)');
    end
    
    % update console
    if opts.verbose && ~opts.dbg
        iterBar.print(ii, delta);
    end

    % update annealing
    if da.curIter < opts.annealingIter
        da.curIter = da.curIter + 1;
    else
        da.pointer = min(da.pointer + 1, da.n);
        da.cur     = da.invbetas(da.pointer);
        da.curIter = 1;
    end
    
    
end

if ~converged && opts.verbose
    iterBar.finish;
    fprintf('EM did not converge in %d iterations\n', opts.maxiter);
end

llh = llh(2:ii+1);
end


function A = stabiliseA_constraintGeneration(obj, curEstimate, verbose)
    % -- what about x0 estimates?
    m          = obj.infer.smooth.mu;
%     P          = obj.infer.smooth.sigma;
%     G          = obj.infer.smooth.G;
%     sumP       = sum(cat(3,P{1:obj.d.T}),3);
%     sumC       = 0;
%     for tt = 1:obj.d.T-1
%         sumC = sumC + P{tt+1} * G{tt}';
%     end    
% 
%     sumC      = sumC;
    s         = obj.suffStats(struct('bIgnoreHash', true));
    S1        = m(:,1:end-1);
    S2        = m(:,2:end);
    
    metric    = obj.par.Q; % doesn't seem to make a difference...?
    metric    = eye(obj.d.x);
%     sumP      = (sumP + sumP)./2;
    A         = ds.utils.learnCGModel_EMQ(S1, S2, s.PHI, s.C', metric, curEstimate, verbose);
end
