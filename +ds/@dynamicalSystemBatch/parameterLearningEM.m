function [llh, ii] = parameterLearningEM(obj, opts)

if nargin < 2 || isempty(opts); opts = struct; end
optsDefault     = struct('epsilon', 1e-3, 'maxiter', 200, 'ssid', false, 'ssidL', 5, ...
                        'verbose', true, 'dbg', false, 'validation', false, ...
                        'multistep', 4, 'diagA', false, 'diagQ', false, 'diagR', false, ...
                        'diagAconstraints', [-1, 1], 'fixBias2', false, ...
                        'sampleStability', 1, 'stableVerbose', false, ...
                        'annealingSchedule', Inf, 'annealingIter', 10, ...
                        'annealingMin', 1e-6, 'strictNegativeCheck', false, ...
                        'filterType', 'linear', 'utpar', struct, 'fixX0', true);
optsDefault     = utils.base.parse_argumentlist(obj.opts, optsDefault, false);      % bring in global opts
opts            = utils.base.parse_argumentlist(opts, optsDefault, false);          % add user specified opts.

sing_val_eps    = 0.004999999999;  % tolerance of singular values of A > 1. This value
                                   % is given in in Sidiqqi et als code for constraint
                                   % generation.
DSS             = obj.allNonEmptyDS;

% ________ Set initial values ____________________________________________
if opts.ssid
    if opts.verbose
        fprintf('(%s) Initialising using subspace identification...\n', datestr(now, 'HH:MM:SS'));
    end
    obj.ssid(opts.ssidL);
    
    obj.smooth('','','forceFilter', true);
%     % do filter and smooth first anyway, since need suff stats for constrain gen
%     if max(abs(eig(obj.par.A))) > 1 + sing_val_eps
%         obj.par.A    = stabiliseA_constraintGeneration(obj, obj.par.A, 1);
%         obj.filter;
%         obj.smooth;
%     end
else
    for dsObj = DSS'
        var_y   = var(dsObj.y);
        if isempty(dsObj.par.A)
            dsObj.par.A = eye(dsObj.d.x);
        end
        if isempty(dsObj.par.Q)
            dsObj.par.Q = var_y*eye(dsObj.d.x)/10;
        end
        if isempty(dsObj.par.H)
            dsObj.par.H = eye(dsObj.d.x);
        end
        if isempty(dsObj.par.R)
            dsObj.par.R = eye(dsObj.d.y);
        end
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
fOpts      = struct('bDoValidation', false, 'bIgnoreHash', true, 'forceFilter', true, 'doLlh', true);
mstepOpts  = struct('verbose', opts.dbg, 'diagQ', opts.diagQ, 'diagR', opts.diagR, ...
                 'diagA', opts.diagA, 'diagAconstraints', opts.diagAconstraints, 'fixBias2', opts.fixBias2);
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
llh        = [-Inf; NaN(opts.maxiter,1)];
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

if opts.diagA; obj.par.A = diag(diag(obj.par.A)); end
if opts.diagQ; obj.par.Q = diag(diag(obj.par.Q)); end
if opts.diagR; obj.par.R = diag(diag(obj.par.R)); end

% Initialise dbg struct in case we need it..
dbgLLH = struct('A',[0,0],'Q',[0,0],'H',[0,0],'R',[0,0],'x0',[0, -Inf]);
bestpar    = cell(1,3);   bestpar{2} = -Inf;

%% MAIN EM LOOP
for ii = 1:opts.maxiter
    % E-Step!
    obj.smooth(opts.filterType, opts.utpar, fOpts);
    
    % llh calc
    llh(ii+1) = obj.logLikelihood('useExisting', true);
    dbgLLH.x0 = [llh(ii+1) - dbgLLH.R(2), llh(ii+1)];
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
        negativeLlhStep = delta < -1e-8 && true && ~(prevStepWasConstrained && ~(opts.sampleStability == 1));
        %negativeLlhStep = delta < -1e-8 && da.cur == 1 && ~(prevStepWasConstrained && ~(opts.sampleStability == 1));
    end
    
    if negativeLlhStep && ii > 1
        iterBar.updateText([iterBar.text, '*']);
        if multiStep == 1
            % basically things have gone really wrong by here..
            iterBar.clearConsole;
%             fprintf('min eigv Q = %.3e', min(eig(obj.par.Q)));
%             fprintf(', min eigv R = %.3e\n', min(eig(obj.par.R)));
            keyboard
        elseif opts.sampleStability > 1
            % only periodically sampling stability of A leads to jumps..
            opts.sampleStability = 1;
        else
            % Doing ECM is not *guaranteed* to increase the llh on each step
            multiStep = multiStep/2;
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
    
    % ---------------------------------------------------------------------
    % ____ Canonical parameters: A ________________________________________
    % cannot do diagonal A and unconstrained B in M-step as yet.
    prevA         = obj.par.A;
    if ~opts.diagA && multiStep > 1
        obj.parameterLearningMStep({'A', 'B'}, mstepOpts);
    else
        % randomly update either A or B first
        for ss = randperm(2)
            if ss == 1;                      obj.parameterLearningMStep({'A'}, mstepOpts); end
            if ss == 2 && obj.hasControl(1); obj.parameterLearningMStep({'B'}, mstepOpts); end
            if multiStep == 1 && obj.hasControl(1); obj.smooth(opts.filterType, opts.utpar, fOpts); end
        end
    end
    
    % Check for stability of A: significant problems when A blows up.
    prevStepWasConstrained = false;
    if mod(ii, opts.sampleStability) == 0 && ~opts.diagA
        error('constrained (non-diagonal) A not implemented yet for dynamicalSystemBatch. Soz.');
    end
    
    if multiStep == 1
        obj.smooth(opts.filterType, opts.utpar, fOpts);
        % if A was constrained, B may no longer be optimal => do over
        if prevStepWasConstrained
            obj.parameterLearningMStep({'B'}, mstepOpts);
            obj.smooth(opts.filterType, opts.utpar, fOpts);
        end
    end
    cLlh      = obj.logLikelihood('useExisting', true);
    dbgLLH.A  = [cLlh - dbgLLH.x0(2), cLlh];
    % ____ (END: parameter A) _____________________________________________
    

    % ____ Canonical parameters: Q, H, R (linear) _________________________
    if multiStep == 4 && obj.emiLinear && obj.evoLinear
        obj.parameterLearningMStep({'Q','H','R'}, mstepOpts);
    elseif obj.evoLinear
        obj.parameterLearningMStep({'Q'}, mstepOpts);
        obj.smooth(opts.filterType, opts.utpar, fOpts);
        cLlh      = obj.logLikelihood('useExisting', true);
        dbgLLH.Q  = [cLlh - dbgLLH.A(2), cLlh];

        obj.parameterLearningMStep({'H'}, mstepOpts);
        if multiStep < 4        
            obj.smooth(opts.filterType, opts.utpar, fOpts);
        end
        cLlh      = obj.logLikelihood('useExisting', true);
        dbgLLH.H  = [cLlh - dbgLLH.Q(2), cLlh];
    
        obj.parameterLearningMStep({'R'}, mstepOpts);
        % likelihood calculation happens later
    else
        % 
    end
    % ____ (END: parameter Q, H, R) _______________________________________
    
    % M-step: x0
    if ~opts.fixX0
        if multiStep == 1
            obj.smooth(opts.filterType, opts.utpar, fOpts); 
            cLlh         = obj.logLikelihood('useExisting', true);
            dbgLLH.R     = [cLlh - dbgLLH.H(2), cLlh];
        else
            dbgLLH.R     = [0, dbgLLH.H(2)];
        end
        
        obj.parameterLearningMStep({'x0'}, mstepOpts);
    else
        dbgLLH.R     = [0, dbgLLH.H(2)];
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
    
    if obj.logLikelihood < bestpar{2}
        obj.par = bestpar{1};
        if opts.verbose
            fprintf('Best log likelihood found on iteration %d of %d\n', bestpar{3}, ii);
        end
    end
        
    obj.smooth(opts.filterType, opts.utpar, fOpts);
end

llh = llh(2:ii+1);
end



function [A, optCode] = stabiliseA_constraintGeneration(obj, prvEstimate, verbose)
    s         = obj.suffStats(struct('bIgnoreHash', true));
    m         = obj.infer.smooth.mu;
    S1        = m(:,1:end-1);
    S2        = m(:,2:end);
    
    metric    = obj.par.Q; % doesn't seem to make a difference...?
%     metric    = eye(obj.d.x);
    [A, optCode]  = ds.utils.learnCGModel_EMQ(S1, S2, s.PHI, s.C', metric, obj.par.A, prvEstimate, verbose);
    if isequal(A, prvEstimate) && optCode >= 0
        optCode = 9999;  % stable matrix unchanged
    end
        
end

function msg = stabiliseA_errorMsg(code)
    if code == 9999
        msg = 'No stable matrix found better than previous estimate.';
        return
    end
    msg       = utils.optim.optimMessage(code, 'onlyErrors', true);
end