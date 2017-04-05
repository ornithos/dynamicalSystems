function [llh, ii] = parameterLearningEM(obj, opts)

if nargin < 2 || isempty(opts); opts = struct; end
optsDefault     = struct('epsilon', 1e-3, 'maxiter', 200, 'ssid', false, 'ssidL', 5, ...
                        'verbose', true, 'dbg', false, 'validation', false, ...
                        'multistep', 4, 'diagA', false, 'diagQ', false, 'diagR', false, ...
                        'diagAconstraints', [-1, 1], 'fixBias2', false, ...
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
% get options for (1) Filtering and (2) M-step
fOpts      = struct('bDoValidation', false, 'bIgnoreHash', true);
mstepOpts  = struct('verbose', opts.dbg, 'diagQ', opts.diagQ, 'diagR', opts.diagR, ...
                 'diagA', opts.diagA, 'diagAconstraints', opts.diagAconstraints, 'fixBias2', opts.fixBias2);
optFds     = fieldnames(opts);
for oname = optFds   % fixA, fixQ, fix....
    if strcmp(oname(1:3), 'fix')
        mstepOpts.(oname) = opts.(oname);
    end
end

% Optimisation options
optimDisplay       = 'none';    % 'iter-detailed'
% optimType          = 'analytic';  % 'auto', 'analytic', 'debug-analytic'

optimOpts          = optimoptions('fminunc','Algorithm','quasi-newton','Display', optimDisplay);
switch opts.optimType
    case 'analytic'
        optimOpts = optimoptions(optimOpts, 'SpecifyObjectiveGradient',true); 
    case 'analytic-debug'
        optimOpts = optimoptions(optimOpts, 'CheckGradients',true, 'FiniteDifferenceType', 'central');
    case 'auto'
        optimOpts = optimoptions(optimOpts, 'FiniteDifferenceType', 'forward');
    otherwise
        error('unknown optimType specified. Choose from ''auto'', ''analytic'', ''analytic-debug''');
end

eljOpts            = struct('gamma', opts.gamma);
optimEmi.options   = optimOpts;
optimEmi.objective = @(x) ds.utilIONLDS.derivEmiWrapper(obj, x, eljOpts);
optimEmi.solver    = 'fminunc';
optimEmi.x0        = [obj.par.emiNLParams.eta(:); obj.par.emiNLParams.C(:)];
LARGE_OPTIM_ITER   = 400;
SMALL_OPTIM_ITER   = 200;

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

% if checking that LLH increases after every M-step part, need to remove
% shortcuts -- these remove guarantees of strictly non-decreasing
if opts.strictNegativeCheck
    opts.multistep = 1;
    opts.sampleStability = 1;
end

% _______ initialise _____________________________________________________
llh        = [-Inf; NaN(opts.maxiter,1)];
converged  = false;
iterBar    = utils.base.objIterationBar;
iterBar.newIterationBar('EM Iteration: ', opts.maxiter, true, '--- ', 'LLH change: ');


% multistep: do all 4 (default) maximisations together initially.
multiStep  = opts.multistep;   % shorter version for readability
% Initialise dbg struct in case we need it..
dbgLLH     = struct('A',[0,0],'Q',[0,0],'H',[0,0],'R',[0,0],'x0',[0,-Inf]);
bestpar    = cell(1,3);

% PARAMETER INTIALISATION: (ensure llh does not drop on ii = 2)
%
% initialise A to be stable if not already. Reqd for constraint generation
% to make meaningful comparison to prevA.
errConstrA = 0;
if (~isfield(mstepOpts,'fixA') || ~mstepOpts.fixA) && opts.sampleStability < opts.maxiter && ...
        max(abs(eig(obj.par.A))) > 1 + sing_val_eps
    obj.par.A     = stabiliseA_constraintGeneration(obj, obj.par.A, zeros(obj.d.x), 0);
    prevA         = obj.par.A;
else
    prevA         = obj.par.A;
end

if opts.diagA; obj.par.A = diag(diag(obj.par.A)); end
if opts.diagQ; obj.par.Q = diag(diag(obj.par.Q)); end
if opts.diagR; obj.par.R = diag(diag(obj.par.R)); end

guaranteedStable = mstepOpts.diagA && max(abs(mstepOpts.diagAconstraints)) <= 1;

%% MAIN EM LOOP
for ii = 1:opts.maxiter
    % E-Step!
    obj.filter('ukf', true, [], fOpts);
    obj.smooth('ukf', [], fOpts);
    
    % llh calc
    dbgLLH.x0  = [obj.infer.llh - dbgLLH.R(2), obj.infer.llh];
    llh(ii+1) = obj.infer.llh;
    delta     = llh(ii+1) - llh(ii);
    
    % ============== Convergence and Admin ===============================
    % negative update warning
    % --- if prev step constrained, we may have seen a drop in LLH.
    %     However, we do not expect this in the case where sampleStability
    %     is 1 since should still be monotonic if every iter stable.
    
    % Llh might increase if constraint generation not used each time, and
    % when it is, it has to revert to previous constrained A. Remove.
    if (errConstrA == 9999) && (opts.sampleStability > 1)
        delta       = delta - dbgLLH.A(1);
        dbgLLH.A(1) = 0;
    end
    
    if ii > 2 && delta < -20
        keyboard
    end
    
    if opts.strictNegativeCheck
        % ensure monotone ON EVERY parameter maximisation
        monotoneLlhFail = dbgLLH.R(1) < -1e-8 || dbgLLH.Q(1) < -1e-8 || dbgLLH.A(1) < -1e-1 || dbgLLH.H(1) < 1e-8;
    else
        % more generous usual monotone checking
        monotoneLlhFail = delta < -1e-8 && ~(prevStepWasConstrained && ~(opts.sampleStability == 1));
    end
    
    if monotoneLlhFail && ii > 1
        iterBar.updateText([iterBar.text, '*']);
        if multiStep == 1
            % basically things have gone really wrong by here..
            iterBar.clearConsole;
            fprintf('min eigv Q = %.3e', min(eig(obj.par.Q)));
            fprintf(', min eigv R = %.3e\n', min(eig(obj.par.R)));
%             keyboard
%             if opts.dbg; keyboard; end
        elseif opts.sampleStability > 1
            % only periodically sampling stability of A leads to jumps..
            opts.sampleStability = 1;
        else
            % Doing ECM is not *guaranteed* to increase the llh on each step
            multiStep = multiStep ./ 2;
        end
    else
        % keep track of best llh found so far (no guarantees of monotonicity
        if llh(ii+1) > bestpar{2}
            bestpar{1} = obj.par;
            bestpar{2} = llh(ii+1);
            bestpar{3} = ii;
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
    % ====================================================================
    %
    %% M-Steps (interleaved with E-steps if required, since ow. ECM not EM!)
    % ----------------------------
    
    % ____ Canonical parameters ___________________________________________
    
%     prevLLH       = obj.logLikelihood;
    obj.parameterLearningMStep({'A', 'B'}, mstepOpts);
%     obj.parameterLearningMStep({'B'}, mstepOpts);
    
    % Check for stability of A: significant problems when A blows up.
    prevStepWasConstrained = false;
    if mod(ii, opts.sampleStability) == 0 && ~guaranteedStable
        if max(abs(eig(obj.par.A))) > 1 + sing_val_eps
            cgTic        = tic;
            if opts.verbose
                iterBar.clearConsole;
                fprintf('Stabilising A with constraint generation... ');
            end
            badA         = obj.par.A;
            [obj.par.A, errConstrA]  = stabiliseA_constraintGeneration(obj, badA, prevA, opts.stableVerbose);
            err          = stabiliseA_errorMsg(errConstrA);
            if opts.verbose
                fprintf('Done! (%.2f)\n', toc(cgTic)); 
                if ~isempty(err); fprintf('STABLE OPTIM: %s\n', err); end
            end
            prevStepWasConstrained = true;
        end
        prevA         = obj.par.A;
    else
        errConstrA    = 0;
    end
    
        if multiStep == 1       
            obj.filter('ukf', true, [], fOpts);
            obj.smooth('ukf', [], fOpts);
            dbgLLH.A  = [obj.infer.llh - dbgLLH.x0(2), obj.infer.llh];
            if prevStepWasConstrained
                obj.parameterLearningMStep({'B'}, mstepOpts);
                obj.filter('ukf', true, [], fOpts);
                obj.smooth('ukf', [], fOpts);
                dbgLLH.A  = [obj.infer.llh - dbgLLH.x0(2), obj.infer.llh];
            end
        end
%     obj.parameterLearningMStep({'B'}, mstepOpts);
%     if obj.logLikelihood < dbgLLH.A(2)
%         iterBar.clearConsole;
%         fprintf('B matrix reduced llh by %.4f\n', obj.logLikelihood - dbgLLH.A(2));
%     end
    
    obj.parameterLearningMStep({'Q'}, mstepOpts);
        if multiStep < 4        
            obj.filter('ukf', true, [], fOpts);
            obj.smooth('ukf', [], fOpts);
            if multiStep==1; dbgLLH.Q  = [obj.infer.llh - dbgLLH.A(2), obj.infer.llh]; end
            if dbgLLH.Q(1) < 0
                if opts.dbg; keyboard; end
            end
        end
       
    % ========= Optimise Nonlinear emission function ===================
    if ~strcmp(optimDisplay, 'none')
        fprintf('\n');
        iterBar.currOutputLen = 0;
    end
    
      if ii <= 3
          optimOpts = optimoptions(optimOpts, 'Display', 'off', 'MaxFunEvals', LARGE_OPTIM_ITER, 'GradObj', 'off');
      elseif ii < opts.maxiter
          optimOpts = optimoptions(optimOpts, 'Display', 'off', 'MaxFunEvals', SMALL_OPTIM_ITER, 'GradObj', 'off');
      else
          optimOpts = optimoptions(optimOpts, 'Display', 'off', 'MaxFunEvals', LARGE_OPTIM_ITER, 'GradObj', 'off');
      end
      optimEmi.options = optimOpts;
      
    if opts.dbg; [F,~,q] = obj.expLogJoint('freeEnergy', true); end
    emiOptOut     = fminunc(optimEmi);  % <- magic happens here
    optimEmi.x0   = emiOptOut;
    obj.par.emiNLParams.eta   = reshape(emiOptOut(1:obj.d.y*4), obj.d.y, 4);
    obj.par.emiNLParams.C     = reshape(emiOptOut((obj.d.y*4+1):end), obj.d.y, obj.d.x);
    if opts.dbg
        [F1,~,q1] = obj.expLogJoint('freeEnergy', true);
        fprintf('M-Step: ''H'' --- Change in FreeNRG: %5.8f\n', F1 - F);
    end

        if multiStep < 4  
            obj.filter('ukf', true, [], fOpts);
            obj.smooth('ukf', [], fOpts);
            if multiStep==1; dbgLLH.H  = [obj.infer.llh - dbgLLH.Q(2), obj.infer.llh]; end
            if dbgLLH.H(1) < -1e-3
                SMALL_OPTIM_ITER = min(SMALL_OPTIM_ITER + 50, LARGE_OPTIM_ITER);
                if opts.verbose
                    iterBar.clearConsole;
                    fprintf('Increased function evals for BFGS + 50\n');
                end
            end
        end
    
    % ========= Calculate R ============================================
    if opts.dbg; [F,~,q] = obj.expLogJoint; end
    [~,M2]        = ds.utilIONLDS.utTransform_ymHx(obj);
    obj.par.R     = M2 ./ obj.d.T;
    if opts.diagR; obj.par.R = diag(diag(obj.par.R));  end
    if multiStep==1
        obj.filter('ukf', true, [], fOpts);
        obj.smooth('ukf', [], fOpts);
        dbgLLH.R  = [obj.infer.llh - dbgLLH.H(2), obj.infer.llh]; 
    end
    if opts.dbg
        [F1,~,q1] = obj.expLogJoint;
        fprintf('M-Step: ''R'' --- Change in FreeNRG: %5.8f\n', F1 - F);
    end
    % ========= (:?) INITIAL DISTN =====================================
    % Not entirely comfortable about this
%     obj.par.x0.mu    = obj.infer.smooth.x0.mu;
%     obj.par.x0.sigma = obj.infer.smooth.x0.sigma;
    
    
    % update console
    if opts.verbose && ~opts.dbg
        iterBar.print(ii, delta);
    end
    lastpar   = obj.par;
end

if ~converged
    if opts.verbose
        iterBar.finish;
        fprintf('EM did not converge in %d iterations\n', opts.maxiter);
    end
    
    if obj.logLikelihood < bestpar{2}
        obj.par = bestpar{1};
        if opts.verbose
            fprintf('Best log likelihood found on iteration %d of %d\n', bestpar{3}, ii);
        end
    end
        
    obj.filter('ukf', true, [], fOpts);
    obj.smooth('ukf', [], fOpts); 
    
    % llh calc
    dbgLLH.x0  = [obj.infer.llh - dbgLLH.R(2), obj.infer.llh];
    llh(ii+1) = obj.infer.llh;
    delta     = llh(ii+1) - llh(ii);
end

llh = llh(2:ii+1);
end


function [A, optCode] = stabiliseA_constraintGeneration(obj, curEstimate, prevEstimate, verbose)
    s         = obj.suffStats(struct('bIgnoreHash', true));
    m         = obj.infer.smooth.mu;
    S1        = m(:,1:end-1);
    S2        = m(:,2:end);
    
    metric    = obj.par.Q; % doesn't seem to make a difference...?
%     metric    = eye(obj.d.x);
    [A, optCode]  = ds.utils.learnCGModel_EMQ(S1, S2, s.PHI, s.C', metric, curEstimate, prevEstimate, verbose);
    if isequal(A, prevEstimate) && optCode >= 0
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