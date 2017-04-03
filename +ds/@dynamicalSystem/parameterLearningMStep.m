function parameterLearningMStep(obj, updateOnly, opts)
% parameterLearningMStep(obj)
% M-Step of the EM algorithm: maximise likelihood of parameters given the
% current estimate of the latent state x.
% This follows the algorithm in Särkkä.
% The updates are exact and analytic: no iteration is required.
%
% Dual updates for [A, B] and [H, C] are constructed by appending the
% control input to x. This is trivially derived, but may be found (eg. in
% Sidiqqi's constraint gen paper (2007)).
%
% The formulation used here is for
%    x_t = Ax_{t-1} + Bu_t.
% Note that it is the control input for the *current* step, not the
% previous one. This is more natural, but both variants can be found, and
% in my notes, the derivations are actually for the control inputs in the
% previous step. This has been altered by incrementing the u's in the
% pairwise expectations.

    if nargin < 2
        updateOnly = [];
    end

    if nargin < 3 || isempty(opts)
        opts = struct;
    end

    if isempty(updateOnly)
        updateOnly = {'A','B','Q','H','C','R'};
    else
        updateOnly = upper(updateOnly);
    end
    
    % obtain fixed quantities, and remove from update list.
    origUpdates = updateOnly;
    optFds      = fieldnames(opts);
    for oo = 1:numel(optFds)   % fixA, fixQ, fix....
        oname = optFds{oo};
        if strcmp(oname(1:3), 'fix')
            rmUpdate = upper(oname(4:end));
            % must search for validity seperately to remove statement,
            % since M-step generally does not include all possible vars.
            if ~ismember(rmUpdate, {'A','B','C','H','Q','R'})
                warning('Unknown fix command: ''fix%s'' in M-step options', rmUpdate);
            end
            findel   = ismember(updateOnly, rmUpdate);
            if opts.(oname) && any(findel)
                updateOnly(findel) = [];
            end
            opts = rmfield(opts, oname);
        end
    end
    
    if isempty(updateOnly)
        if opts.verbose
            fprintf('M-step: No non-fixed quantities in %s to update. Exiting...\n', strjoin(origUpdates, ', '));
        end
        return
    end
    
    optsDefault = struct('verbose', true, 'diagQ', false, 'diagR', false, ...
                         'diagA', false, 'diagAconstraints', [-1,1], 'anneal', 1);
    opts        = utils.base.parse_argumentlist(opts, optsDefault, true);
    
    % get all sufficient statistics
    s   = obj.suffStats(struct('bIgnoreHash', true, 'anneal', opts.anneal));

    eps = 0;% 1e-8./obj.d.T;   % 'kind-of' prior to avoid instability in covariances.
    
    % opts to stop log likelihood doing unnecessary work
    fOpts = struct('bDoValidation', false, 'bIgnoreHash', true);

%% Parameter updates
    if opts.verbose; prevLLH = obj.logLikelihood('',[],fOpts); end
    
    % ______________ A (And B) updates ___________________________________
    matchesAB = ismember({'A','B'}, updateOnly);
    if any(matchesAB)
        % A / B updates where A is unconstrained - a full matrix
        if ~opts.diagA
            if obj.hasControl(1)
                if all(matchesAB)
                    AB          = [s.C, s.XU] / [s.PHI, s.Xm1_U; s.Xm1_U', s.UU];
                    obj.par.A   = AB(:, 1:obj.d.x);
                    obj.par.B   = AB(:, obj.d.x+1:end);
                elseif matchesAB(1)
                    obj.par.A   = (s.C -  obj.par.B * s.Xm1_U') / s.PHI;
                elseif matchesAB(2)
                    obj.par.B   = (s.XU - obj.par.A * s.Xm1_U) / s.UU;
                end
            elseif matchesAB(1)
                obj.par.A   = s.C / s.PHI;
            end
        else
        % A / B updates where A is constrained to be diagonal. This is not
        % the same as projecting optimum of A onto the space of diag mats.
            if matchesAB(1)
                psi = inv(obj.par.Q);   % .. not ideal
                if ~obj.hasControl(1)
                    tmp = s.C';
                else
                    tmp = s.C' -  s.Xm1_U * obj.par.B';
                end
                v   = diag(tmp * psi);  %#ok (<- preproc) numericall less stable, but more efficient.
                M   = s.PHI .* psi;
                % solve
                a   = M\v;
                
                % ? possible projection constraints
                if numel(opts.diagAconstraints) == 2
                    con = sort(opts.diagAconstraints);
                    % if violating, must solve as linearly constrained least squares
                    if any(a>con(2)) || any(a<con(1))
                        lb            = repmat(con(1), obj.d.x, 1); ub = repmat(con(2), obj.d.x, 1);
                        a             = min(ub, max(a, lb));
                        % ------------ EXACT MINIMISATION -----------------
%                         lsqopts       = optimoptions('fmincon', 'Display', 'none');
%                         minfn         = @(x) -2*trace(psi*diag(x)*tmp) + trace(psi*diag(x)*s.PHI*diag(x));
%                         [aOpt,~,xfl]  = fmincon(minfn,a,[],[],[],[],lb,ub,[],lsqopts);
%                         if xfl < 0; a = diag(obj.par.A); end
                        % -------------------------------------------------
                        % Greedy active-set heuristic -- much faster.. haven't encountered any problems
                        % yet, but may fail occasionally. --> if worried use the exact version above.
                        afix = NaN(size(a));
                        while numel(a) > 0 && (any(a < lb) || any(a > ub))
                            cur                = find(isnan(afix));
                            aViolLB            = a < lb;
                            afix(cur(aViolLB)) = lb(aViolLB);
                            v                  = v - M(:,cur(aViolLB)) * lb(aViolLB);
                            aViolUB            = a > ub;
                            afix(cur(aViolUB)) = ub(aViolUB);
                            v                  = v - M(:,cur(aViolUB)) * ub(aViolUB);
                            M(:,cur(aViolLB | aViolUB)) = [];
                            a                  = M \ v;
                        end
                        afix(isnan(afix))  = a;
                        a                  = afix;
%                         fprintf('normdiff aOpt to aActive = %.7f\n', norm(aOpt - a));
                    end
                elseif ~isempty(opts.diagAconstraints)
                    error(isempty(opts.diagAconstraints), 'diagAconstraints must be a 2 element vector or empty');
                end
                obj.par.A = diag(a);
            end
            if matchesAB(2) && obj.hasControl(1)
                obj.par.B   = (s.XU - obj.par.A * s.Xm1_U) / s.UU;
            end
        end
    end
    
    % debug message
    if opts.verbose; prevLLH = dbgConsole(prevLLH, obj, updateOnly(ismember(updateOnly, {'A','B'})), fOpts); end
    
    
    % ______________ Q updates ___________________________________________
    if any(strcmpi(updateOnly, 'Q'))
        % original Q update (without inputs)
        A           = obj.par.A;
        Q           = s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A';
        
        if obj.hasControl(1)
            % --- additions to covariance from control inputs ----
            B       = obj.par.B;
            Bum     = B * s.XU';
            BuAm_m  = B * s.Xm1_U' * A';
            Q       = Q + B*s.UU*B' - Bum - Bum' + BuAm_m + BuAm_m';
        end

        % numerical imprecision (hopefully!) on symmetry 
        obj.par.Q   = (Q + Q')./2 + eps*eye(obj.d.x);
        
        if opts.diagQ
            obj.par.Q = diag(diag(obj.par.Q));
        end
        
        % debug message
        if opts.verbose; prevLLH = dbgConsole(prevLLH, obj, {'Q'}, fOpts); end
    end

    % ______________ H (And C) updates ___________________________________
    matchesHC = ismember({'H','C'}, updateOnly);
    if any(matchesHC)
        s2 = s.emissions;
        
        HCopts = struct;
        HCopts.H       = matchesHC(1);
        HCopts.C       = matchesHC(2) && obj.hasControl(2);  % clearly only do C update if control exists.
        HCopts.bias    = ~isempty(obj.par.c);
        HCopts.control = obj.hasControl(2);
        
        % find "numerator" and "denominator" and calculate regression
        if HCopts.H
            if HCopts.C
                if HCopts.bias
                    % H && C && bias 
                    numer     = [s.B, s.YU, s2.Ymu];
                    denom     = [s2.SIGMA, s2.XU, s2.Xmu;   s2.XU', s2.UU, s2.Umu;   s2.Xmu', s2.Umu', 1];
                    result    = numer / denom;
                    obj.par.H = result(:, 1:obj.d.x);
                    obj.par.C = result(:, obj.d.x+1:end-1);
                    obj.par.c = result(:, end);
                else
                    % H && C && ____
                    numer     = [s.B, s.YU];
                    denom     = [s2.SIGMA, s2.XU;   s2.XU', s2.UU];
                    result    = numer / denom;
                    obj.par.H = result(:, 1:obj.d.x);
                    obj.par.C = result(:, obj.d.x+1:end-1);
                end
            else
                if HCopts.bias
                    if HCopts.control
                        % H && _ && bias && control
                        numer     = [(s.B - obj.par.C * s2.XU'), s2.Ymu];
                        denom     = [s2.SIGMA, s2.Xmu;  s2.Xmu', 1];
                        result    = numer / denom;
                        obj.par.H = result(:, 1:obj.d.x);
                        obj.par.c = result(:, end);
                    else
                        % H && _ && bias && ______
                        numer     = [s.B, s2.Ymu];
                        denom     = [s2.SIGMA, s2.Xmu;  s2.Xmu', 1];
                        result    = numer / denom;
                        obj.par.H = result(:, 1:obj.d.x);
                        obj.par.c = result(:, end);
                    end
                else
                    if HCopts.control
                        % H && _ && ____ && control
                        numer     = s.B - obj.par.C * s2.XU';
                        denom     = s2.SIGMA;
                        result    = numer / denom;
                        obj.par.H = result(:, 1:obj.d.x);
                    else
                        % H && _ && ____ && ______
                        numer     = s.B;
                        denom     = s2.SIGMA;
                        result    = numer / denom;
                        obj.par.H = result(:, 1:obj.d.x);
                    end
                end
            end
        else
                if HCopts.bias
                    % H && C && bias 
                    numer     = [(s.YU - obj.par.H * s2.XU), s2.Ymu];
                    denom     = [s2.UU, s2.Umu;   s2.Umu', 1];
                    result    = numer / denom;
                    obj.par.C = result(:, 1:obj.d.u);
                    obj.par.c = result(:, end);
                else
                    % H && C && ____
                    numer     = s.YU - obj.par.H * s2.XU;
                    denom     = s2.UU;
                    result    = numer / denom;
                    obj.par.C = result(:, 1:obj.d.u);
                end
        end
        
    end
    
    % debug message
    if opts.verbose; prevLLH = dbgConsole(prevLLH, obj, updateOnly(ismember(updateOnly, {'H','C'})), fOpts); end
    
    % ______________ R updates ___________________________________________
    if any(strcmpi(updateOnly, 'R'))
        % original R update (without inputs)
        s2          = s.emissions;
        H           = obj.par.H;
        R           = s.D - H*s.B' - s.B*H' + H*s2.SIGMA*H';
        
        if obj.hasControl(2)
            if any(any(isnan(obj.y))); warning('NaNs not tested properly in emission learning'); end
            % --- additions to covariance from control inputs ----
            C       = obj.par.C;
            Cuy     = C * s2.YU';
            Cum_H   = C * s2.XU' * H';
            R       = R + C*s2.UU*C' - Cuy - Cuy' + Cum_H + Cum_H';
            if ~isempty(obj.par.c)
                % bias exists && control exists.
                tmp = obj.par.c*s2.Umu'*obj.par.C';
                R   = R + (-tmp - tmp')./s2.Ty;
            end
        end

        if ~isempty(obj.par.c)
            % bias exists
            tmp = obj.par.c*(s2.Ymu - obj.par.H * s2.Xmu);
            R   = R + (-tmp - tmp' + obj.par.c*obj.par.c')./s2.Ty;
        end
        
        % numerical imprecision (hopefully!) on symmetry
        obj.par.R   = (R + R')./2 + eps*eye(obj.d.y);
        
        if opts.diagR
            obj.par.R = diag(diag(obj.par.R));
        end
        
        % debug message
        if opts.verbose; prevLLH = dbgConsole(prevLLH, obj, {'R'}, fOpts); end
    end
 
end


function newLLH = dbgConsole(prevLLH, obj, label, fOpts) 
    newLLH = obj.logLikelihood('',[],fOpts); 
    fprintf('M-Step: ''%s'' --- Change in LLH: %5.8f\n', strjoin(label, ''','''), newLLH - prevLLH);
end