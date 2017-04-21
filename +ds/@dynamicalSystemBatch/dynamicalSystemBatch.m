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
                    yyHat(:,tt)    = obj.doEmission(xx(:,tt+1), u_t);
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
                        cFit(:,tt) = obj.doEmission(obj.infer.smooth.mu{nn}(:,tt), u_t);
                    else
                        cFit(:,tt) = obj.doEmission(obj.infer.smooth.mu{nn}(:,tt), u_t, obj.infer.smooth.sigma{nn}{tt}, utpar);
                    end
                end
                fitted{nn} = cFit;
            end
            if nargout == 0
                % overwrite existing values.
                obj.yhat = fitted;
            end
        end

        function fitted = getPredictedValues(obj, nlookahead, utpar)
            if nargin < 2; nlookahead = 0; end
            obj.ensureInference('PREDVALS', 'filter');
        
            if nlookahead == -Inf
                fitted = getFittedValues(obj);
                return;
            end
            assert(nlookahead >= 0, 'lagged smoothing values not implemented. Try nlookahead=-Inf for Smoothed');
            
            %if ~obj.emiLinear; [~,~,h,Dh]     = obj.functionInterfaces; nlpars.f = h; nlpars.Df = Dh; end
            doLinOrEKF  = nargin < 3 || isempty(utpar);
            
            fitted = cell(obj.d.n, 1);
            for nn = 1:obj.d.n
                cFit      = zeros(obj.d.y, obj.d.T(nn) - nlookahead);
                for tt = 1:obj.d.T - nlookahead
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
                        cFit(:,tt) = obj.doEmission(x_t, u_t);
                    else
                        cFit(:,tt) = obj.doEmission(x_t, u_t, P, utpar);
                    end
                end
                fitted{nn} = cFit;
            end
        end
        
        function predict = getPredictFreeRun(obj, t, l, utpar)
            % predict from (time t), (l datapoints) forward
            if nargin < 2; t = 1; end
            assert(utils.is.scalarint(t) && t > 0 && t <= obj.d.T, 't must be a scalar int in 1,...,T');
            if nargin < 3
            	if t == obj.d.T; error('length of output l must be specified'); end
                l = obj.d.T - t; 
            end
            assert(utils.is.scalarint(l) && l > 0, 'l must be a positive scalar int');
            obj.ensureInference('PREDVALS', 'filter');

            doLinOrEKF = nargin < 4 || isempty(utpar);
            
            predict = cell(obj.d.n, 1);
            for nn = 1:obj.d.n
                if ~doLinOrEKF; P = obj.infer.filter.sigma{nn}{t+1}; end

                cPred = NaN(obj.d.y, t+l);
                x_t = obj.infer.filter.mu{nn}(:, t);
                for tt = t+1:t+l
                    if tt <= obj.d.T && any(obj.hasControl)
                        u_t = obj.u(:,tt); 
                    else
                        u_t = 0;
                    end
                    x_t                 = obj.doTransition(x_t, u_t);
                    if ~doLinOrEKF; P = obj.par.A*P*obj.par.A' + obj.par.Q; end
                    if doLinOrEKF
                        cPred(:,tt) = obj.doEmission(x_t, u_t);
                    else
                        cPred(:,tt) = obj.doEmission(x_t, u_t, P, utpar);
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
end