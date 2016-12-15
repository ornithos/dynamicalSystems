function [llh, ii] = parameterLearningEM(obj, opts)

if nargin < 2 || isempty(opts); opts = struct; end
optsDefault     = struct('epsilon', 1e-3, 'maxiter', 200, 'ssid', false, 'ssidL', 5, ...
                        'verbose', true, 'dbg', false, 'validation', false, ...
                        'multistep', 4, 'diagQ', false, 'diagR', false, ...
                        'sampleStability', 1, 'stableVerbose', false, 'optimType', 'analytic');
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


% ________ Get names of parameters which have been 'fixed' _______________
% get options for (1) Filtering and (2) M-step
fOpts      = struct('bDoValidation', false, 'bIgnoreHash', true);
mstepOpts  = struct('verbose', opts.dbg, 'diagQ', opts.diagQ, 'diagR', opts.diagR);
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
optimEmi.options   = optimOpts;
optimEmi.objective = @(x) ds.utilIONLDS.derivEmiWrapper(obj, x);
optimEmi.solver    = 'fminunc';
optimEmi.x0        = [obj.par.emiNLParams.eta(:); obj.par.emiNLParams.C(:)];

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

% _______ initialise _____________________________________________________
llh        = [-Inf; zeros(opts.maxiter,1)];
converged  = false;
iterBar    = utils.base.objIterationBar;
iterBar.newIterationBar('EM Iteration: ', opts.maxiter, true, '--- ', 'LLH change: ');


% multistep: do all 4 (default) maximisations together initially.
multiStep  = opts.multistep;
% Initialise dbg struct in case we need it..
dbgLLH = struct('A',[0,0],'Q',[0,0],'H',[0,0],'R',[0,-Inf]);

%% MAIN EM LOOP
for ii = 1:opts.maxiter
    % E-Step!
    obj.filter('ukf', true, [], fOpts);
    obj.smooth('ukf', [], fOpts);
    
    % llh calc
    dbgLLH.R  = [obj.infer.llh - dbgLLH.H(2), obj.infer.llh];
    llh(ii+1) = obj.infer.llh;
    delta     = llh(ii+1) - llh(ii);
    
    % ============== Convergence and Admin ===============================
    % negative update warning
    % --- if prev step constrained, we may have seen a drop in LLH.
    %     However, we do not expect this in the case where sampleStability
    %     is 1 since should still be monotonic if every iter stable.
    if delta < -1e-8 && ~(prevStepWasConstrained && ~(opts.sampleStability == 1))
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
    % ====================================================================
    %
    %% M-Steps (interleaved with E-steps if required, since ow. ECM not EM!)
    % ----------------------------
    
    % ____ Canonical parameters ___________________________________________
    prevA         = obj.par.A;
%     prevLLH       = obj.logLikelihood;
    obj.parameterLearningMStep({'A','B'}, mstepOpts);
%     obj.parameterLearningMStep({'B'}, mstepOpts);
    
    % Check for stability of A: significant problems when A blows up.
    prevStepWasConstrained = false;
    if mod(ii, opts.sampleStability) == 0
        if max(abs(eig(obj.par.A))) > 1 + sing_val_eps
            cgTic        = tic;
            if opts.verbose
                iterBar.clearConsole;
                fprintf('Stabilising A with constraint generation... ');
            end
            badA         = obj.par.A;
            obj.par.A    = prevA;   % reset to last good A
            obj.par.A    = stabiliseA_constraintGeneration(obj, prevA, opts.stableVerbose); % see mini function below ????
%             newLLH       = obj.logLikelihood;
%             if newLLH < prevLLH
%                 keyboard
%             end
            if opts.verbose; fprintf('Done! (%.2f)\n', toc(cgTic)); end
            prevStepWasConstrained = true;
        end
    end
    
        if multiStep == 1       
            obj.filter('ukf', true, [], fOpts);
            obj.smooth('ukf', [], fOpts);
            dbgLLH.A  = [obj.infer.llh - dbgLLH.R(2), obj.infer.llh];
        end
    
    obj.parameterLearningMStep({'Q'}, mstepOpts);
        if multiStep < 4        
            obj.filter('ukf', true, [], fOpts);
            obj.smooth('ukf', [], fOpts);
            if multiStep==1; dbgLLH.Q  = [obj.infer.llh - dbgLLH.A(2), obj.infer.llh]; end
        end
       
    % ========= Optimise Nonlinear emission function ===================
    if ~strcmp(optimDisplay, 'none')
        fprintf('\n');
        iterBar.currOutputLen = 0;
    end
    
      optimOpts = optimoptions(optimOpts, 'MaxFunEvals',30,'GradObj', 'off');
      if ii <= 3
          optimOpts = optimoptions(optimOpts, 'Display', 'iter-detailed', 'MaxFunEvals', 200, 'GradObj', 'off');
      end
      if ii == opts.maxiter
          optimOpts = optimoptions(optimOpts, 'MaxFunEvals', 100, 'GradObj', 'off');
      end
      optimEmi.options = optimOpts;
      
    emiOptOut     = fminunc(optimEmi);  % <- magic happens here
    optimEmi.x0   = emiOptOut;
    obj.par.emiNLParams.eta   = reshape(emiOptOut(1:obj.d.y*5), obj.d.y, 5);
    obj.par.emiNLParams.C     = reshape(emiOptOut((obj.d.y*5+1):end), obj.d.y, obj.d.x);
    
        if multiStep < 4  
            obj.filter('ukf', true, [], fOpts);
            obj.smooth('ukf', [], fOpts);
            if multiStep==1; dbgLLH.H  = [obj.infer.llh - dbgLLH.Q(2), obj.infer.llh]; end
        end
    
    % ========= Calculate R ============================================
    [~,M2]        = ds.utilIONLDS.utTransform_ymHx(obj);
    obj.par.R     = M2 ./ obj.d.T;
    if opts.diagR; obj.par.R = diag(diag(obj.par.R)); end
    
    % ========= (:?) INITIAL DISTN =====================================
    % Not entirely comfortable about this
    obj.par.x0.mu    = obj.infer.smooth.x0.mu;
    obj.par.x0.sigma = obj.infer.smooth.x0.sigma;
    
    
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
