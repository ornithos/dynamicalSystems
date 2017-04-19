function [llh, ii] = parameterLearningEM_bspl(obj, opts)

if nargin < 2 || isempty(opts); opts = struct; end
optsDefault     = struct('epsilon', 1e-3, 'maxiter', 200, 'ssid', false, 'ssidL', 5, ...
                        'verbose', true, 'dbg', false, 'validation', false, ...
                        'multistep', 4, 'diagA', false, 'diagQ', false, 'diagR', false, ...
                        'diagAconstraints', [-1, 1], 'fixBias2', false, ...
                        'sampleStability', 1, 'stableVerbose', false, ...
                        'annealingSchedule', Inf, 'annealingIter', 10, ...
                        'annealingMin', 1e-6, 'strictNegativeCheck', false, ...
                        'filterType', 'linear', 'utpar', struct, 'fixX0', true, ...
                        'nonlinOptimOpts', struct);
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
fOpts      = struct('bDoValidation', false, 'bIgnoreHash', true, 'forceFilter', true, 'doLlh', true);
mstepOpts  = struct('verbose', opts.dbg, 'diagQ', opts.diagQ, 'diagR', opts.diagR, ...
                 'diagA', opts.diagA, 'diagAconstraints', opts.diagAconstraints, 'fixBias2', opts.fixBias2);
optFds     = fieldnames(opts);
for oname = optFds   % fixA, fixQ, fix....
    if strcmp(oname(1:3), 'fix')
        mstepOpts.(oname) = opts.(oname);
    end
end

%% _________ Optimisation options ________________________________________

nonlinOpts                  = struct;
nonlinOpts.optimDisplay     = 'none';   %'iter-detailed'; %'none';
nonlinOpts.large_iter       = 400;      % initial # BFGS iterations
nonlinOpts.small_iter       = 100;      % subsequent # iterations - only small mvmts required.
nonlinOpts.optimType        = 'fminunc';
nonlinOpts.chgEtas          = true(1, numel(obj.par.emiNLParams.bspl.t));  % default change all etas
nonlinOpts.etaUB            = [];
nonlinOpts.bfgsSpline       = false;
nonlinOpts.monotoneWarning  = true;
nonlinOpts                  = utils.struct.structCoalesce(opts.nonlinOptimOpts, nonlinOpts);


% transform eta --> theta (hat). Used for optimisation object x0.
theta              = [obj.par.emiNLParams.eta(:,1), diff(obj.par.emiNLParams.eta, [], 2)];
if any(theta(:) <= 0)
    if nonlinOpts.monotoneWarning; warning('some values of eta are nonincreasing. Replacing where relevant with diff of 1e-12'); end
    theta = max(theta, 1e-12);
end
u                  = log(theta);
adaptU             = u(:, nonlinOpts.chgEtas);


% Set up optimisation object for fminunc / fmincon.
assert(ismember(nonlinOpts.optimType, {'fminunc', 'fmincon'}), 'optimType must be ''fminunc'' or ''fmincon''');
optimOpts          = optimoptions(nonlinOpts.optimType, 'Display', nonlinOpts.optimDisplay, ...
                                    'SpecifyObjectiveGradient', true, 'CheckGradients', false, 'MaxFunEvals', nonlinOpts.large_iter);
if strcmp(nonlinOpts.optimType, 'fminunc')
    optimOpts = optimoptions(optimOpts, 'Algorithm','quasi-newton');
end
optimEmi           = struct;
optimEmi.options   = optimOpts;
optimEmi.solver    = nonlinOpts.optimType;
optimEmi.x0        = [adaptU(:); obj.par.emiNLParams.C(:); obj.par.emiNLParams.bias(:)];
% (----> OPTIMISATION OBJECTIVE <-----)
optimEmi.objective = @(x) ds.utilIONLDS.derivEmiWrapper_bspl(obj, x, 'fixBias2', opts.fixBias2, ...
                        'etaMask', nonlinOpts.chgEtas, 'etaUB', nonlinOpts.etaUB, 'bfgsSpline', nonlinOpts.bfgsSpline);

% Setup constraints for spline coefficients
nonlinOpts.A       = kron(eye(obj.d.y), ones(1,sum(nonlinOpts.chgEtas)));
nonlinOpts.b       = nonlinOpts.etaUB;
if strcmp(nonlinOpts.optimType, 'fmincon')
    optimEmi.Aineq         = [nonlinOpts.A, zeros(obj.d.y, obj.d.y*(obj.d.x + 1))];
    optimEmi.bineq         = nonlinOpts.b;
end

% for ease of reference, add chgEtas to workspace.
chgEtas            = nonlinOpts.chgEtas;

%%
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
    dbgLLH.x0  = [obj.infer.llh - dbgLLH.R(2), obj.infer.llh];
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
        if max(abs(eig(obj.par.A))) > 1 + sing_val_eps
            cgTic        = tic;
            if opts.stableVerbose > 1
                iterBar.clearConsole;
                fprintf('Stabilising A with constraint generation... ');
            end
%             badA         = obj.par.A;
%             obj.par.A    = prevA;   % reset to last good A
            [obj.par.A, errConstrA]  = stabiliseA_constraintGeneration(obj, prevA, opts.stableVerbose); % see mini function below 
            err          = stabiliseA_errorMsg(errConstrA);
            if opts.verbose
                fprintf('Done! (%.2f)\n', toc(cgTic)); 
                if ~isempty(err); fprintf('STABLE OPTIM: %s\n', err); end
            end
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
    
    if multiStep == 1
        obj.smooth(opts.filterType, opts.utpar, fOpts);
        % if A was constrained, B may no longer be optimal => do over
        if prevStepWasConstrained
            obj.parameterLearningMStep({'B'}, mstepOpts);
            obj.smooth(opts.filterType, opts.utpar, fOpts);
        end
    end
    dbgLLH.A  = [obj.infer.llh - dbgLLH.x0(2), obj.infer.llh];
    % ____ (END: parameter A) _____________________________________________
    

    % ____ Canonical parameters: Q ________________________________________
    obj.parameterLearningMStep({'Q'}, mstepOpts);
    obj.smooth(opts.filterType, opts.utpar, fOpts);
    dbgLLH.Q  = [obj.infer.llh - dbgLLH.A(2), obj.infer.llh];
    % ____ (END: parameter Q) _____________________________________________
    
    %% ========= Optimise Nonlinear emission function ===================
    if ~strcmp(nonlinOpts.optimDisplay, 'none')
        fprintf('\n');
        iterBar.currOutputLen = 0;
    end
    
    % set number of optimisation iterations
    if ii <= 3 || ii == opts.maxiter
        optimOpts = optimoptions(optimOpts, 'MaxFunEvals', nonlinOpts.large_iter);
    elseif ii == 4   % reset to small iter
        optimOpts = optimoptions(optimOpts, 'MaxFunEvals', nonlinOpts.small_iter);
    end
    optimEmi.options = optimOpts;
      
    if opts.dbg; [F,~,q] = obj.expLogJoint_bspl('freeEnergy', true); end
    
    % Perform optimisation
    if strcmp(nonlinOpts.optimType , 'fminunc')
        emiOptOut     = fminunc(optimEmi);  % <- magic happens here
    else
        emiOptOut     = fmincon(optimEmi);
    end
    optimEmi.x0   = emiOptOut;
    
    % update parameters from optimisation search
    ds.utilIONLDS.updateParams_bspl(obj, emiOptOut, chgEtas);

%     % test gradients (orig file vs monotonic file are same)
%     [fMono,gradMono] = ds.utilIONLDS.bsplineGradMono(obj, []);
%     [fOrig,gradOrig] = ds.utilIONLDS.bsplineGrad(obj, []);
%     assert(abs(fMono - fOrig) < 1e-9, 'fval differences > 1e-9');
%     assert(max(max(abs(gradMono.eta - gradOrig.eta))) < 1e-10, 'eta grad differences > 1e-10');
%     assert(max(max(abs(gradMono.C - gradOrig.C))) < 1e-10, 'C grad differences > 1e-10');
%     assert(max(max(abs(gradMono.bias - gradOrig.bias))) < 1e-10, 'bias grad differences > 1e-10');
    
    % DEBUG VISUALISATION
%     figure
%     xtmp = bsxfun(@plus, obj.par.emiNLParams.C*obj.infer.smooth.mu, obj.par.emiNLParams.bias);
%     xtmp = floor(min(min(xtmp))):ceil(max(max(xtmp)));
%     for zzz = 1:2
%         subplot(1,2,zzz); plot(xtmp, obj.par.emiNLParams.bspl.functionEval(xtmp, obj.par.emiNLParams.eta(zzz,:))); hold on; 
%         plot(bsxfun(@plus, obj.par.emiNLParams.C(zzz,:)*obj.infer.smooth.mu, obj.par.emiNLParams.bias(zzz)), obj.y(zzz,:), '*'); hold off; 
%     end
    
    % if optimising spline outside of BFGS, solve small QP.
    if ~nonlinOpts.bfgsSpline
        obj.smooth(opts.filterType, opts.utpar, fOpts);
        qpopts                              = optimoptions('quadprog', 'Display', 'none');
        eta0                                = exp(reshape(emiOptOut(1:sum(chgEtas)*obj.d.y),obj.d.y,sum(chgEtas)))';
        [~,~,mm]                            = ds.utilIONLDS.bsplineGradMono(obj, emiOptOut, opts.utpar, 'etaMask', chgEtas);
        [eta,~,xfl]                         = quadprog(mm.K, -mm.v, nonlinOpts.A, nonlinOpts.b, [],[],...
                                                        zeros(size(eta0(:))),[],eta0(:), qpopts);
        utils.optim.optimMessage(xfl, 'onlyErrors', true);
        obj.par.emiNLParams.eta(:,chgEtas)  = cumsum(reshape(eta, sum(chgEtas), obj.d.y), 1)';
        tfeta                               = reshape(log(eta), sum(chgEtas), obj.d.y)';
        optimEmi.x0(1:obj.d.y*sum(chgEtas)) = tfeta(:);
    
    
%         figure
%         for zzz = 1:2
%             subplot(1,2,zzz); plot(xtmp, obj.par.emiNLParams.bspl.functionEval(xtmp, obj.par.emiNLParams.eta(zzz,:))); hold on; 
%             plot(bsxfun(@plus, obj.par.emiNLParams.C(zzz,:)*obj.infer.smooth.mu, obj.par.emiNLParams.bias(zzz)), obj.y(zzz,:), '*'); hold off; 
%         end
    end
    
    if opts.dbg
        [F1,~,q1] = obj.expLogJoint_bspl('freeEnergy', true);
        fprintf('M-Step: ''H'' --- Change in FreeNRG: %5.8f\n', F1 - F);
    end
    
        if multiStep < 4  
            obj.smooth(opts.filterType, opts.utpar, fOpts);
            if multiStep==1; dbgLLH.H  = [obj.infer.llh - dbgLLH.Q(2), obj.infer.llh]; end
            if dbgLLH.H(1) < -1e-3
                nonlinOpts.small_iter = min(nonlinOpts.small_iter + 50, nonlinOpts.large_iter);
                if opts.verbose
                    iterBar.clearConsole;
                    fprintf('Increased function evals for BFGS + 50\n');
                end
            end
        end
        
    % ____ Canonical parameters: R ________________________________________
    if opts.dbg; [F,~,q] = obj.expLogJoint; end
    [~,M2]        = ds.utilIONLDS.utTransform_ymHx_bspl(obj);
    
    % NOT SUPER HAPPY ABOUT THIS: SHOULD I REMOVE ALL COMPLETELY MISSING Ys
    % FROM M2, and DIVIDE BY Ty?
    obj.par.R     = M2 ./ obj.d.T;                     
    if opts.diagR; obj.par.R = diag(diag(obj.par.R));  end
    if multiStep==1
        obj.smooth(opts.filterType, opts.utpar, fOpts);
        dbgLLH.R  = [obj.infer.llh - dbgLLH.H(2), obj.infer.llh]; 
    end
    if opts.dbg
        [F1,~,q1] = obj.expLogJoint;
        fprintf('M-Step: ''R'' --- Change in FreeNRG: %5.8f\n', F1 - F);
    end
    
    % M-step: x0
    if ~opts.fixX0
        if multiStep == 1
            obj.smooth(opts.filterType, opts.utpar, fOpts); 
            dbgLLH.R     = [obj.infer.llh - dbgLLH.H(2), obj.infer.llh];
        else
            dbgLLH.R     = [0, dbgLLH.H(2)];
        end
        
        obj.par.x0.mu    = obj.infer.smooth.x0.mu;
        obj.par.x0.sigma = obj.infer.smooth.x0.sigma;
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