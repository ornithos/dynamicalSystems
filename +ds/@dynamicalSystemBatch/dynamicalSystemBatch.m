classdef dynamicalSystemBatch < ds.dynamicalSystem
    % dynamicalSystemBatch: create batch version of dynamicalSystems object
    %
    % dynamicalSystemBatch(dimx, dimy, 'evolution', {A || f || []}, {Df, Q}
    %                 'emission', {H || h || []}, {Dh, R},
    %                 'data', {y || T}, 'x0' {x0mu, x0cov}, 
    %                 ('xtrue', x),
    %                 ('control', u, (C), (D)),
    %                 opts)
    %
    % dynamicalSystemBatch({dsObj1, dsObj2, ...});
    %
    % Each argument is the same as for dynamicalSystem, but is optionally a
    % cell containing multiple copies of time series with assumed same
    % parameters.
    
    properties (Hidden = true)
        sglObj = []
    end
    
    properties (Dependent)
        intrinsicDimension
    end
    
    methods
        function obj = dynamicalSystemBatch(varargin)
            
            % hack as need to do superclass constructor first..
%             dsVarargin = varargin;
%             for ii = 1:numel(dsVarargin)
%                 if iscell(dsVarargin{ii})
%                     dsVarargin{ii} = dsVarargin{ii}{1};
%                 end
%             end
            obj@ds.dynamicalSystem(varargin{:});
            
            % handle copy!
            if numel(varargin) == 1 && isa(varargin{1}, 'ds.dynamicalSystemBatch')
                return
            end
        end

        
        function yy = generateData(obj)
            contyn       = utils.base.questionUser('This will overwrite existing data. Continue', {'y','n'}, 'caseSensitive', false);
            if strcmpi(contyn, 'n'); return; end
            
            cellX = cell(obj.d.n,1); cellY = cell(obj.d.n, 1); cellYH = cell(obj.d.n,1);
            transChol    = chol(obj.par.Q)';
            emissChol    = chol(obj.par.R)';
            for nn = 1:obj.d.n
                xx       = zeros(obj.d.x, obj.d.T(nn)+1);
                yy       = zeros(obj.d.y, obj.d.T(nn));
                xx(:,1)  = obj.par.x0.mu;
                yyHat    = zeros(obj.d.y, obj.d.T(nn));
            
                for tt = 1:obj.d.T(nn)
                    u_t = [];
                    if any(obj.hasControl); u_t = obj.u{nn}(:,tt); end
                    xx(:,tt+1)     = obj.doTransition(xx(:,tt), u_t) + transChol * randn(obj.d.x,1);
                    yyHat(:,tt)    = obj.doEmission(xx(:,tt+1), u_t, [], [], nn);
                    yy(:,tt)       = yyHat(:,tt) + emissChol * randn(obj.d.y,1);
                end
                
                cellX{nn} = xx(:,2:end); cellY{nn} = yy; cellYH{nn} = yyHat;
            end
             
            % overwrite existing values.
            obj.x        = cellX;
            obj.y        = cellY;
            obj.yhat     = cellYH;
        end
        
        % Make a copy of a handle object.
        function new = copy(this)
            % Instantiate new object of the same class.
            curWarns           = this.opts.warnings;
            this.opts.warnings = false;
            new                = ds.dynamicalSystemBatch(this);
            this.opts.warnings = curWarns;
        end
        
        function fitted = getFittedValues(obj, utpar)
            obj.ensureInference('FITVALS', 'smooth');
            fitted = cell(obj.d.n, 1);
            for nn = 1:obj.d.n
                cFit = zeros(obj.d.y, obj.d.T(nn));
                for tt = 1:obj.d.T(nn)
                    u_t = [];
                    if any(obj.hasControl); u_t = obj.u{nn}(:,tt); end
                    if nargin < 2 || isempty(utpar)
                        cFit(:,tt) = obj.doEmission(obj.infer.smooth.mu{nn}(:,tt), u_t, [], [], nn);
                    else
                        cFit(:,tt) = obj.doEmission(obj.infer.smooth.mu{nn}(:,tt), u_t, obj.infer.smooth.sigma{nn}{tt}, utpar, nn);
                    end
                end
                fitted{nn} = cFit;
            end
            if nargout == 0
                % overwrite existing values.
                obj.yhat = fitted;
            end
        end

        function [fitted, X] = getPredictedValues(obj, nlookahead, utpar, nRng)
            if nargin < 2; nlookahead = 0; end
            obj.ensureInference('PREDVALS', 'filter');
        
            if isinf(nlookahead)
              if nlookahead < 0; fitted = obj.getFittedValues; end
              if nlookahead > 0; fitted = obj.getPredictFreeRun(1); end
              return;
            end
            
            assert(nlookahead >= 0, 'lagged smoothing values not implemented. Try nlookahead=-Inf for Smoothed');
            
            %if ~obj.emiLinear; [~,~,h,Dh]     = obj.functionInterfaces; nlpars.f = h; nlpars.Df = Dh; end
            doLinOrEKF  = nargin < 3 || isempty(utpar);
            
            if nargin < 4 || isempty(nRng); nRng = 1:obj.d.n; end
            fitted     = cell(obj.d.n, 1); 
            X          = cell(obj.d.n, 1);
            
            for nn = nRng
                cFit      = zeros(obj.d.y, obj.d.T(nn) - nlookahead);
                cX        = zeros(obj.d.x, obj.d.T(nn) - nlookahead);
                for tt = 1:obj.d.T(nn) - nlookahead
                    x_t = obj.infer.filter.mu{nn}(:,tt);
                    u_t = [];
                    if any(obj.hasControl); u_t = obj.u{nn}(:,tt); end
                    if ~doLinOrEKF; P = obj.infer.filter.sigma{nn}{tt}; end
                    for jj = 1:nlookahead
                        x_t = obj.doTransition(x_t, u_t);
                        u_t = [];
                        if any(obj.hasControl); u_t = obj.u{nn}(:,tt+jj); end
                        if ~doLinOrEKF; P = obj.par.A*P*obj.par.A' + obj.par.Q; end
                    end
                    if doLinOrEKF
                        cFit(:,tt) = obj.doEmission(x_t, u_t, [], [], nn);
                    else
                        cFit(:,tt) = obj.doEmission(x_t, u_t, P, utpar, nn);
                    end
                    cX(:,tt) = x_t;
                end
                fitted{nn} = cFit;
                X{nn}      = cX;
            end
        end
        
        function predict = getPredictFreeRun(obj, t, l, utpar)
            % predict from (time t), (l datapoints) forward
            if nargin < 2; t = 1; end
            assert(all(arrayfun(@(x) utils.is.scalarint(x, 1), t)) && all(t <= obj.d.T), 't must be scalar in 1, ..., T');
            if numel(t) == 1 && obj.d.n > 1; t = repmat(t, obj.d.n, 1); end
            if nargin < 3
            	if all(t == obj.d.T); error('length of output l must be specified'); end
                l = obj.d.T - t; 
            else
                if numel(l) == 1; l = repmat(l, obj.d.n, 1); end
                assert(numel(l) == obj.d.n, 'l must have either 1 or %d elements', obj.d.n);
            end
            assert(all(arrayfun(@(x) utils.is.scalarint(x, 1), l)), 'l must be a positive scalar int');
            obj.ensureInference('PREDVALS', 'filter');

            doLinOrEKF = nargin < 4 || isempty(utpar);
            
            predict = cell(obj.d.n, 1);
            for nn = 1:obj.d.n
                if ~doLinOrEKF; P = obj.infer.filter.sigma{nn}{t(nn)+1}; end

                cPred = NaN(obj.d.y, t(nn)+l(nn));
                x_t = obj.infer.filter.mu{nn}(:, t(nn));
                for tt = t(nn)+1:t(nn)+l(nn)
                    if tt <= obj.d.T(nn) && any(obj.hasControl)
                        u_t = obj.u{nn}(:,tt); 
                    else
                        u_t = 0;
                    end
                    x_t                 = obj.doTransition(x_t, u_t);
                    if ~doLinOrEKF; P = obj.par.A*P*obj.par.A' + obj.par.Q; end
                    if doLinOrEKF
                        cPred(:,tt) = obj.doEmission(x_t, u_t, [], [], nn);
                    else
                        cPred(:,tt) = obj.doEmission(x_t, u_t, P, utpar, nn);
                    end
                end
                predict{nn} = cPred;
            end
        end
        
        % not yet implemented stuff -- overwrite so don't call ds fns.
        function removeControl(obj)   %#ok (args not used - prototype)
            error('control manipulation functions not implemented in batch.');
        end
        function addControl(obj, u, evo, emi) %#ok (args not used - prototype)
            error('control manipulation functions not implemented in batch.');
        end
        function truncateTime(obj, T) %#ok (args not used - prototype)
            error('time manipulation functions not implemented in batch.');
        end
        
        %% DYNAMICAL SYSTEM METHODS       
        function ssid(obj, varargin)   %#ok (args not used - prototype)
            error(['SSID not implemented for multiple time series yet. ', ....
                'See Holcomb & Bitmead 2017 (https://arxiv.org/pdf/1704.02635.pdf)']);
        end
        
        
        
        % prototype
%         [llh, niters] = parameterLearningEM(obj, opts);
        D             = filter(obj, fType, bDoLLH, utpar, opts);
                        smooth(obj, sType, utpar, varargin)
        s             = suffStats(obj, varargin);
        [y, covY]     = impute_y(obj, varargin);
        
        % ***********************************************************
        % --> NEED TO IMPLEMENT PARAMETERLEARNINGMSTEP (INCL x0!!) FOR BTCH
        % ***********************************************************
      
        
        %% MISC FUNCTIONS
        
        function value = get.intrinsicDimension(obj)
            value = cellfun(@(x) sum(~all(isnan(x),2)), obj.y);
        end
        function set.intrinsicDimension(obj, val)   %#ok
            % do nothing! Critical for on copy properties.
        end
        
        function dsObj = extractSingle(obj, nn, doInfer)
            assert(utils.is.scalarint(nn,1) && nn <= obj.d.n, 'nn is not scalar int in 1, .., N');
            if nargin < 3 || isempty(doInfer); doInfer = true; end
            
            % construct prototype object if not exist yet
            if isempty(obj.sglObj)
                opts       = struct('warnings', false, 'verbose', false);
                if isa(obj, 'ds.ionlds')
                    
                    [~,~,h,~] = obj.functionInterfaces;
                    obj.sglObj = ds.ionlds(obj.d.x, obj.d.y, 'x0', obj.par.x0.mu{1}, ...
                        obj.par.x0.sigma{1}, 'evolution', eye(obj.d.x), eye(obj.d.x), 'emission', ...
                        h, eye(obj.d.y), 'data', obj.y{1}, 'control', obj.u{1}, obj.par.B, false, opts);
                    obj.sglObj.par.emiNLParams = obj.par.emiNLParams;
                    
                elseif isa(obj, 'ds.dynamicalSystem')
                    obj.sglObj = ds.dynamicalSystem(obj.d.x, obj.d.y, 'x0', obj.par.x0.mu{1}, ...
                        obj.par.x0.sigma{1}, 'evolution', eye(obj.d.x), eye(obj.d.x), 'emission', ...
                        ones(obj.d.y, obj.d.x), eye(obj.d.y), 'data', obj.y{1}, 'control', obj.u{1}, false, false, opts);
                end
                for i = obj.sglObj.stackptr:-1:1; obj.sglObj.delete(i); end
                obj.sglObj.infer = ds.internal.emptyInferStruct;
            end

            % Copy all properties across: move from batch --> single by
            % extracting nn'th element from each cell (except stack and
            % xxxx) which are not applicable in this case).
            p      = properties(obj.sglObj);
            dsObj  = obj.sglObj.copy;
            % valid while properties of dsBatch are a superset of ds:
            for i = 1:length(p)
                if ismember(p{i}, {'stack'}); continue; end   % infer
                val          = obj.(p{i});
                cur          = dsObj.(p{i});
                if isstruct(cur)
                    dsObj.(p{i}) = ds.dynamicalSystemBatch.copyStructNoCellRecursive(dsObj.(p{i}), val, nn);
                else
                    if iscell(val); val = val{nn}; end
                    dsObj.(p{i}) = val;
                end
            end
            
            % unfortunately, obj.d.T is a *MATRIX* for batch, and scalar ow
            dsObj.d.T = obj.d.T(nn);
            
            if doInfer
                dsObj.smooth(dsObj.infer.fType, [], 'forceFilter', true);
            end
        end
%         
%         function out = n(obj)
%             out = numel(obj.dsArray);
%         end
% 
%         function idx = firstNonEmptyIdx(obj)
%             idx = find(~cellfun(@isempty, obj.dsArray), 1, 'first');
%         end
% 
%         function arr = allNonEmptyDS(obj)
%             idx = ~cellfun(@isempty, obj.dsArray);
%             arr = obj.dsArray(idx);
%         end
%         
%         function garbageCollect(obj)
%             elEmpty      = cellfun(@isempty, obj.dsArray);
%             obj.dsArray(elEmpty) = [];
%         end
% 
%         function ok  = validationInference(obj, varargin)
%             dsarr  = obj.allNonEmptyDS;
%             ok     = true;
%             for ii = 1:numel(dsarr)
%                 ok = ok && dsarr{ii}.validationInference(varargin{:});
%             end
%         end
% 
%         % index access to object
%         function sref = subsref(obj,s)
%             switch s(1).type
%                 case '.'
%                     noOutput = {'filter', 'smooth', 'ssid', 'garbageCollect'};
%                     if ismember(s(1).subs, noOutput)
%                         builtin('subsref', obj, s);
%                     else
%                         sref = builtin('subsref',obj,s);
%                     end
% %                 case '()'
% %                     sref = builtin('subsref',obj,s);
%                 case '{}'
%                     % in case of multiple subscripts > pass remaining to ds
%                     ll = length(s);
% 
%                     curS = s(1);
%                     curS.type = '{}';
%                     assert(numel(curS.subs) == 1, 'dynamicalSystemBatch object is 1D');
% 
%                     sref = builtin('subsref', obj.dsArray, curS);
% 
%                     if ll > 1
%                         curS = s(2:end);
%                         sref = builtin('subsref', sref, curS);                     
%                     end
% 
%                 otherwise
%                     error('bad suscript type: %s (-- I only accept {})', s(1).type);
%             end
%         end

    end

    methods (Access = public, Hidden=true)
        parameterLearningMStep(obj, updateOnly, opts); % internals for EM
        
        % ---- Dynamics wrappers --------------------------
%         function out = doTransition(obj, input, u)
%            
%             if obj.evoLinear
%                   out        = obj.par.A * input;
%                   if obj.hasControl(1); out = out + obj.par.B * u; end
%             else
%                 [f,~,~,~]    = obj.functionInterfaces;
%                 if ~obj.hasControl(1), out = f(input);
%                 else, out = f(input, u); end
%             end
%         end
        
        function out = doEmission(obj, input, u, P, utpar, nn)
            
            % get bias (may be null, vector, or cell)
            if ~isempty(obj.par.c)
                if iscell(obj.par.c) && nargin == 6 && ~isempty(nn)
                    b = obj.par.c{nn}; 
                else
                    b = obj.par.c;
                end
            else
                b  = zeros(obj.d.y, 1);
            end
            
            if obj.emiLinear
                 out = obj.par.H * input + b;
                 if obj.hasControl(2); out = out + obj.par.C * u; end
            else
                [~,~,h,~]    = obj.functionInterfaces;
                if nargin >= 5 && ~isempty(utpar)
                    nlpars.f = h;
                    nlpars.Q = 0;
                    out   = ds.utils.assumedDensityTform(nlpars, input, P, u, 2, utpar);
                else
                    if ~obj.hasControl(2); out   = h(input);
                    else, out = h(input, u); end
                end
            end
        end
        % ------------------------------------------------
    end
   
    methods (Access = protected)

    end
    
    methods (Static)
        function s1 = copyStructNoCellRecursive(s1, s2, nn)
            % copies struct s2 into s1, traversing struct as a tree.
            % --> the nn'th value is taken from all cells at all levels of
            % the struct.
            %
            % ASSUMES ALL FIELDS WITHIN STRUCT s1 ARE AVAILABLE IN s2.
            % (NO CHECKS - WILL BREAK OW.)
            
            if isempty(s1)
                s1 = s2;
                return
            end
            
            fdss1 = fieldnames(s1);
            for ii = 1:numel(fdss1)
                f   = fdss1{ii};
                val = s2.(f);
                if isstruct(val)
                    s1.(f) = ds.dynamicalSystemBatch.copyStructNoCellRecursive(s1.(f), val, nn);
                    continue
                elseif iscell(val)
                    val = val{nn}; 
                end
                s1.(f) = val;
            end
        end
    end
end