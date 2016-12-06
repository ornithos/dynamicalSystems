function parameterLearningMStepNew(obj, updateOnly, opts)
% parameterLearningMStep(obj)
% M-Step of the EM algorithm: maximise likelihood of parameters given the
% current estimate of the latent state x.
% This follows the algorithm in S�rkk�.
% The updates are exact and analytic: no iteration is required.

    if nargin < 2
        updateOnly = [];
    end

    if nargin < 3 || isempty(opts)
        opts = struct;
    end
    optsDefault = struct('verbose', true, 'diagQ', false, 'diagR', false);
    opts        = utils.base.parse_argumentlist(opts, optsDefault, true);

    if isempty(updateOnly)
        updateOnly = {'A','B','Q','H','C','R'};
    else
        updateOnly = upper(updateOnly);
    end
    
    % obtain fixed quantities, and remove from update list.
    origUpdates = updateOnly;
    optFds      = fieldnames(opts);
    for oname = optFds   % fixA, fixQ, fix....
        if strcmp(oname{1}(1:3), 'fix')
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
        end
    end
    
    if isempty(updateOnly)
        if opts.verbose
            fprintf('M-step: No non-fixed quantities in %s to update. Exiting...\n', strjoin(origUpdates, ', '));
        end
        return
    end
            
    % get all sufficient statistics
    s   = obj.suffStats(struct('bIgnoreHash', true));

    eps = 1e-8./obj.d.T;   % 'kind-of' prior to avoid instability in covariances.
    
    % opts to stop log likelihood doing unnecessary work
    fOpts = struct('bDoValidation', false, 'bIgnoreHash', true);

%% Parameter updates
    if opts.verbose; prevLLH = obj.logLikelihood('',[],fOpts); end
    
    % ______________ A (And B) updates ___________________________________
    matchesAB = ismember({'A','B'}, updateOnly);
    if any(matchesAB)
        if obj.hasControl(1)
            if all(matchesAB)
                AB          = [s.C, s.Xm1Um1] / [s.PHI, s.Xm1Um1; s.Xm1Um1', s.Um1Um1];
                obj.par.A   = AB(:, 1:obj.d.x);
                obj.par.B   = AB(:, obj.d.x+1:end);
            elseif matchesAB(1)
                obj.par.A   = (s.C -  obj.par.B * s.Xm1Um1') / s.PHI;
            elseif matchesAB(2)
                obj.par.B   = (s.X_Um1 - obj.par.A * s.Xm1Um1) / s.Xm1Um1;
            end
        elseif matchesAB(1)
            obj.par.A   = s.C / s.PHI;
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
            Bum     = B * s.X_Um1;
            BuAm_m  = B * s.Xm1_U * A';
            Q       = Q + B*s.Um1Um1*B' - Bum - Bum' - BuAm_m - BuAm_m';
        end

        % numerical imprecision on symmetry (hopefully!)
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
        if obj.hasControl(2)
            if all(matchesHC)
                HC          = [s.B, s.XU] / [s.SIGMA, s.YU; s.YU', s.UU];
                obj.par.H   = HC(:, 1:obj.d.x);
                obj.par.C   = HC(:, obj.d.x+1:end);
            elseif matchesHC(1)
                obj.par.H   = (s.B  -  obj.par.C * s.XU') / s.SIGMA;
            elseif matchesHC(2)
                obj.par.C   = (s.YU -  obj.par.H * s.XU)  / s.UU;
            end
        elseif matchesHC(1)
            obj.par.H   = s.B / s.SIGMA;
        end
    end
    
    % debug message
    if opts.verbose; prevLLH = dbgConsole(prevLLH, obj, updateOnly(ismember(updateOnly, {'H','C'})), fOpts); end
    
    % ______________ R updates ___________________________________________
    if any(strcmpi(updateOnly, 'R'))
        % original R update (without inputs)
        H           = obj.par.H;
        R           = s.D - H*s.B' - s.B*H' + H*s.SIGMA*H';
        
        if obj.hasControl(2)
            % --- additions to covariance from control inputs ----
            C       = obj.par.C;
            Cum     = C * s.XU';
            Cum_H   = Cum * H';
            R       = R + C*s.UU*C' - Cum - Cum' - Cum_H - Cum_H';
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