classdef dynamicalSystem < handle
   % dynamicalSystem: create linear/non-linear dynamical system object with
   % access to filtering, smoothing and learning methods, and visualisations
   % for comparison of methods.
   %
   % dynamicalSystem(dimx, dimy, 'evolution', {A || f || []}, {Df, Q}
   %                 'emission', {H || h || []}, {Dh, R},
   %                 'data', {y || T}, 'x0' {x0mu, x0cov}, 
   %                 ('xtrue', x),
   %                 ('control', u, (C), (D)),
   %                 opts)
   %
   % The string arguments can be redistributed as desired, but we assume
   % that they appear in this order for documentation. 
   %
   % dimx - dimension of latent (evolutionary) state.
   % dimy - dimension of output observation (emission).
   %
   % *** EVOLUTION ***
   % - the first argument is either the (linear) transition matrix, a (non-
   % linear) function handle f, or an empty vector. The empty vector is
   % interpreted as an unknown linear matrix which must be learned.
   % - a couple of optional arguments. Df = the derivative (Hessian) of f,
   % the evolution noise covariance Q (can be learned if not specified).
   % - If f is parameterised, the parameters must be entered as a struct
   % into f's final argument. This struct must precede the function handle
   % in the dynamicalSystem inputs.
   %
   % *** EMISSION ***
   % - the first argument is the same as that of the evolution parameters,
   % corresponding to the emission linear/nonlinear function.
   % - a couple of optional arguments. Dh = the derivative (Hessian) of h,
   % the emission noise covariance R (can be learned if not specified).
   % % - If h is parameterised, the parameters must be entered as a struct
   % into h's final argument. This struct must precede the function handle
   % in the dynamicalSystem inputs.
   %
   % *** DATA ***
   % - either a data matrix y corresponding to the observed data (each
   % column is an observation) at time 1, ..., t; or the number of
   % observations that should be generated from the given dynamics. The
   % latter may only be specified if all transitions and covariances have
   % also been specified.
   %
   % ***** xtrue *****  (~~OPTIONAL~~)
   % - x: the true latent state, if we have access to it.
   %
   % ***** x0 *****
   % - x0mu: the mean of the prior over the latent space (if not given
   % assumed 0)
   % - x0cov: the covariance of the prior over the latent space (required).
   %
   % *** control *** (~~OPTIONAL~~)
   % - u: the control inputs (in R^k) as a (T x k) matrix
   % - C: (optional), the control matrix in the latent state
   % - D: (optional), the control matrix in the emission (typically zero).
   % Note that D may be specified without leaving a gap for C if it can be
   % identified from its dimensionality.
   %
   % *** OPTS ***  (~~OPTIONAL~~)
   % - A struct containing fields: warnings, verbose
   % Note that no specifier is given in this case: the final field is 
   % assumed to be opts if it is a struct.
   
   properties
      opts, stack = cell(100,2),
      %  ___ Parameters ___
       par = struct('x0',struct,'A',[],'H',[],'Q',[],'R',[],'c', [], ...
                      'f',[],'Df',[],'h',[],'Dh',[], ...
                      'evoNLParams',struct,'emiNLParams',struct, ...
                      'B',[],'C',[], 'uu', [])
   end
   properties (SetAccess = protected)
       x, y, yhat, d, u
       evoLinear = []
       evoNLhasParams = false
       emiLinear = []
       emiNLhasParams = false
       hasControl = false(2,1);
       stackptr = 0
       % ----------------- Changeable:
       %  ___ Inference ___
       infer = struct('filter',[], 'smooth', [], 'llh', [], 'fType', [], ...
                      'sType', [], 'fpHash', [])
   end
   methods
      function obj = dynamicalSystem(varargin)
         % CONSTRUCTOR
         obj.processInputArgs(varargin);
         
         if numel(varargin) == 1 && (isa(varargin{1}, 'ds.dynamicalSystem') || ...
             strcmp(varargin{1}, 'nullConstruct'))
             return
         end
         
         % pre populate filter/smoother for input data
         if obj.evoLinear && obj.emiLinear
             if obj.validationInference(false)
                 if obj.opts.verbose; fprintf('(%s) Running smoother for input parameters...  ', datestr(now, 'HH:MM:SS')); end
                 obj.filter('Kalman');
                 obj.smooth('Linear');
                 if obj.opts.verbose; fprintf('Complete!\n'); end
                 obj.save('initialised');
             else
                 if obj.opts.warnings
                     warning('dynamicalSystem cannot perform initial inference since not all parameters specified.');
                 end
             end
         else
             if obj.opts.warnings
                 warning('dynamicalSystem has not performed initial inference (nonlinear): choice of algm required.');
             end
         end
         
      end
      
      function yy = generateData(obj)
          xx       = zeros(obj.d.x, obj.d.T+1);
          yy       = zeros(obj.d.y, obj.d.T);
          xx(:,1)  = obj.par.x0.mu;
          yyHat    = zeros(obj.d.y, obj.d.T);
          
          % xtrue given
          if ~isempty(obj.x)
              xx = [xx(:,1), obj.x];
              doX = false;
          else
              doX = true;
          end
          
          transChol   = chol(obj.par.Q)';
          emissChol   = chol(obj.par.R)';
          for tt = 1:obj.d.T
              u_t = [];
              if any(obj.hasControl); u_t = obj.u(:,tt); end
              if doX; xx(:,tt+1)    = obj.doTransition(xx(:,tt), u_t) + transChol * randn(obj.d.x,1); end
              yyHat(:,tt)   = obj.doEmission(xx(:,tt+1), u_t);
              yy(:,tt)      = yyHat(:,tt) + emissChol * randn(obj.d.y,1);
              
              %plot(obj.par.H(1,:) * obj.x(:,1:tt-1), obj.par.H(2,:) * obj.x(:,1:tt-1)); hold on; plot(obj.y(1,1:tt-1), obj.y(2,1:tt-1)); hold off;
              %pause
          end
          if nargout == 0
              % overwrite existing values.
              obj.x       = xx(:,2:end);
              obj.y       = yy;
              obj.yhat    = yyHat;
              obj.infer.fpHash = [];
          end
      end
      
      % Make a copy of a handle object.
      function new = copy(this)
          % Instantiate new object of the same class.
          curWarns           = this.opts.warnings;
          this.opts.warnings = false;
          new                = ds.dynamicalSystem(this);
          this.opts.warnings = curWarns;
      end
      
      function fitted = getFittedValues(obj, utpar)
          obj.ensureInference('FITVALS', 'smooth');
          fitted = zeros(obj.d.y, obj.d.T);
          for tt = 1:obj.d.T
              u_t = [];
              if any(obj.hasControl); u_t = obj.u(:,tt); end
              if nargin < 2 || isempty(utpar)
                fitted(:,tt) = obj.doEmission(obj.infer.smooth.mu(:,tt), u_t);
              else
                fitted(:,tt) = obj.doEmission(obj.infer.smooth.mu(:,tt), u_t, obj.infer.smooth.sigma{tt}, utpar);
              end
          end
          if nargout == 0
              % overwrite existing values.
              obj.yhat = fitted;
          end
      end
      
      function [fitted, X, Cov] = getPredictedValues(obj, nlookahead, utpar, nRng)      %#ok nRng unused; conformity with batch.
          if nargin < 2; nlookahead = 0; end
          if nargin < 3; utpar = []; end
          obj.ensureInference('PREDVALS', 'filter');
          
          if isinf(nlookahead)
              if nlookahead < 0; fitted = obj.getFittedValues(utpar); end
              if nlookahead > 0; [fitted, X] = obj.getPredictFreeRun(1, [], utpar); end
              return;
          end
          assert(nlookahead >= 0, 'lagged smoothing values not implemented. Try nlookahead=-Inf for Smoothed');
          
          if nargout > 2
              saveCovariance = true;
              Cov = cell(obj.d.T - nlookahead,1);
          else
              saveCovariance = false;
          end
          
          %if ~obj.emiLinear; [~,~,h,Dh]    = obj.functionInterfaces; nlpars.f = h; nlpars.Df = Dh; end
          fType      = ds.utils.filterTypeLookup(obj.infer.fType, true) -1;
          parsEmission = ds.internal.getParams(obj, 2, 2*(fType > 0));
          doCov      = fType == 2 || saveCovariance;
          
          fitted     = zeros(obj.d.y, obj.d.T - nlookahead);
          X          = zeros(obj.d.x, obj.d.T - nlookahead);
          
          
          for tt = 1:obj.d.T - nlookahead
              x_t = obj.infer.filter.mu(:,tt);
              u_t = [];
              if any(obj.hasControl); u_t = obj.u(:,tt); end
              if doCov; P = obj.infer.filter.sigma{tt}; end
              for jj = 1:nlookahead
                  x_t = obj.doTransition(x_t, u_t);
                  u_t = [];
                  if any(obj.hasControl); u_t = obj.u(:,tt+jj); end
                  if doCov; P = obj.par.A*P*obj.par.A' + obj.par.Q; end
              end
              if saveCovariance && ~(fType == 2 && isempty(utpar))   % last condition is only for back-compatibility with some old code. Suggest bin it.
                  if ~obj.hasControl(2); u_t = []; end
                  [fm, fP] = ds.utils.assumedDensityTform(parsEmission, x_t, P, [], 2, utpar);
                  fitted(:,tt) = fm;
                  if saveCovariance; Cov{tt} = fP; end
              else
                  fitted(:,tt) = obj.doEmission(x_t, u_t);
              end
              X(:,tt) = x_t;
          end
      end
      
      function [predict, X] = getPredictFreeRun(obj, t, l, utpar)
          % predict from (time t), (l datapoints) forward
          if nargin < 2; t = 1; end
          assert(utils.is.scalarint(t) && t > 0 && t <= obj.d.T, 't must be a scalar int in 1,...,T');
          if nargin < 3 || isempty(l)
              if t == obj.d.T; error('length of output l must be specified'); end
              l = obj.d.T - t; 
          end
          assert(utils.is.scalarint(l) && l > 0, 'l must be a positive scalar int');
          obj.ensureInference('PREDVALS', 'filter');
          
          doLinOrEKF = nargin < 4 || isempty(utpar);
          if ~doLinOrEKF; P = obj.infer.filter.sigma{t+1}; end
          
          predict = NaN(obj.d.y, t+l);
          X       = NaN(obj.d.x, t+1);
          
          x_t = obj.infer.filter.mu(:, t);
          for tt = t+1:t+l
              if tt <= obj.d.T && any(obj.hasControl)
                  u_t = obj.u(:,tt); 
              else
                  u_t = 0;
              end
              x_t             = obj.doTransition(x_t, u_t);
              if ~doLinOrEKF; P = obj.par.A*P*obj.par.A' + obj.par.Q; end
              if doLinOrEKF
                predict(:,tt) = obj.doEmission(x_t, u_t);
              else
                predict(:,tt) = obj.doEmission(x_t, u_t, P, utpar);
              end
              X(:,tt) = x_t;
          end
      end
      
        
      %% ----- Save user interfaces --------------------
      function save(obj, descr)
          assert(nargin == 2, 'please provide a description');
          obj.ensureInference('SAVE', 'smooth');
          
          insertion         = struct;
          insertion.par     = obj.par;
          insertion.infer   = obj.infer;
          insertion.yhat    = obj.yhat;
          insertion.hasControl = obj.hasControl;
          stackPush(obj, insertion, descr);
      end
      
      function delete(obj, descr)
          obj.stackDelete(obj.stackFind(descr));
      end
      
      function descr = savedList(obj)
          sList = obj.getStackDescrList;
          for ii = 1:numel(sList)
              sList{ii} = [num2str(ii), '. ', sList{ii}];
          end
          descr = strjoin(sList, '\n');
      end
      
      function svPoint = getSaved(obj, savedName)
          idx     = obj.stackFind(savedName);
          svPoint = obj.stack{idx, 1};
      end
      
      function useSavedParameters(obj, savedName, verbose)
          if nargin < 3 || isempty(verbose); verbose = obj.opts.verbose; end
          idx     = obj.stackFind(savedName);
          svPoint = obj.stack{idx, 1};
          if verbose; fprintf('Using parameters from save-point ''%s''..\n', obj.stack{idx,2}); end
          obj.par   = svPoint.par;
          obj.infer = svPoint.infer;
          obj.yhat  = svPoint.yhat;
      end
      
      function yhat = getFittedFromSaved(obj, savedName)
          tmpsvName  = char(floor(94*rand(1, 20)) + 32);   % 1e39 possibilities
          obj.save(tmpsvName);
          try
              obj.useSavedParameters(savedName, false);
          catch ME
              obj.stackDelete(obj.stackFind(tmpsvName)); % <- avoid temp saves from failures.
              if strcmp(ME.identifier, 'ds:saveNameNotOnStack')
                  error('save name %s not found. Unable to retrieve fitted values.', savedName);
              else
                  rethrow(ME);
              end
          end
          obj.getFittedValues;
          if isempty(obj.yhat); obj.getFittedValues; end
          yhat    = obj.yhat;
          
          obj.useSavedParameters(tmpsvName, false);
          obj.stackDelete(obj.stackFind(tmpsvName));  % belt-and-braces
      end
      

      %% -----  Other utils ------------------------
      % --- (NL) functions ---------------
      function [f, Df, h, Dh] = functionInterfaces(obj, returnAsStruct)
        % non-linear wrappers for common interface to fns
        if obj.evoNLhasParams
            if obj.hasControl(1)
                f = @(x,u) obj.par.f(x, u, obj.par.evoNLParams);
                Df = @(x,u) obj.par.Df(x, u, obj.par.evoNLParams);
            else
                f = @(x,u) obj.par.f(x, obj.par.evoNLParams);
                Df = @(x,u) obj.par.Df(x, obj.par.evoNLParams);
            end
        else
            if obj.hasControl(1)
                f = @(x,u) obj.par.f(x,u);
                Df = @(x,u) obj.par.Df(x,u);
            else
                f = @(x,u) obj.par.f(x);
                Df = @(x,u) obj.par.Df(x);
            end
        end
        if obj.emiNLhasParams
            if obj.hasControl(2)
                h = @(x,u) obj.par.h(x, u, obj.par.emiNLParams);
                Dh = @(x,u) obj.par.Dh(x, u, obj.par.emiNLParams);
            else
                h = @(x) obj.par.h(x, obj.par.emiNLParams);
                Dh = @(x) obj.par.Dh(x, obj.par.emiNLParams);
            end
        else
            if obj.hasControl(2)
                h = @(x,u) obj.par.h(x,u);
                Dh = @(x,u) obj.par.Dh(x,u);
            else
                h = @(x,u) obj.par.h(x);
                Dh = @(x,u) obj.par.Dh(x);
            end
        end
        if nargin > 1 && returnAsStruct
            assert(nargout <= 1, 'too many outputs requested for returnAsStruct');
            f = struct('f', f, 'Df', Df, 'h', h, 'Dh', Dh);
        end
      end
      
      function b   = parametersChanged(obj)
            b = ~strcmp(obj.infer.fpHash, obj.parameterHash);
      end
      
      function ensureInference(obj, caller, type)
          % ensure inference exists, and if linear, perform inference if not.
          assert(ischar(caller), 'caller must be of type char');
          if nargin < 3; type = 'smooth'; end
          assert(ischar(type) && any(ismember(type, {'filter', 'smooth'})), 'type must be ''filter'' or ''smooth''');
          if obj.parametersChanged
              if obj.opts.warnings; fprintf('%s: posterior does not match current parameters', upper(caller)); end
              if obj.evoLinear && obj.emiLinear
                  if obj.opts.warnings; fprintf('... fixed!\n'); end
                  ftype = obj.infer.fType;
                  stype = obj.infer.sType;
                  obj.filter(ftype);
                  if strcmp(type, 'smooth')
                      obj.smooth(stype);
                  end
              else
                  error('Unable to get fitted values: posterior does not match parameters. Please run a posterior algm');
              end
          end
      end
          
      function llh = logLikelihood(obj, fType, utpar, opts)
          if nargin < 4 || isempty(opts); opts = struct; end
          if nargin < 3 || isempty(utpar); utpar = struct; end
          if nargin < 2 || isempty(fType); fType = obj.infer.sType; end
          tmpobj = obj.copy;
          tmpobj.filter(fType, true, utpar, opts);
          llh    = tmpobj.infer.llh;
      end
      
      function removeControl(obj)
            if ~any(obj.hasControl)
                warning('obj does not have a control to remove');
            end
            obj.hasControl = false(2,1);
            obj.u = [];
      end
      function addControl(obj, u, evo, emi)
            assert(~any(obj.hasControl), 'obj already has control');
            assert(nargin == 4, 'require 3 arguments: u, evo (bool), emi (bool)');
            assert(isscalar(evo) && islogical(evo), 'evo must be scalar boolean');
            assert(isscalar(emi) && islogical(emi), 'emi must be scalar boolean');
            assert(size(u, 2) == obj.d.T, 'u must have T = %d columns', obj.d.T);
            
            obj.hasControl = [evo, emi];
            obj.d.u        = size(u, 2);
            if evo && (isempty(obj.par.B) || all(size(obj.par.B) == [obj.d.x, obj.d.u]))
                obj.par.B      = zeros(obj.d.x, obj.d.u);
            end
            if emi && (isempty(obj.par.C) || all(size(obj.par.C) == [obj.d.y, obj.d.u]))
                obj.par.C      = zeros(obj.d.y, obj.d.u);
            end
            obj.u = u;
      end
      function truncateTime(obj, T)
          assert(utils.is.scalarint(T, 1), 'T must be a positive scalar integer');
          assert(T <= obj.d.T, 'T is greater than the current max. time');
          
          obj.d.T      = T;
          obj.u        = obj.u(:,1:T);
          obj.y        = obj.y(:,1:T);
          if ~isempty(obj.x);    obj.x    = obj.x(:,1:T); end
          if ~isempty(obj.yhat); obj.yhat = obj.yhat(:,1:T); end
          
          obj.infer.filter.mu    = obj.infer.filter.mu(:,1:T);
          obj.infer.filter.sigma = obj.infer.filter.sigma(1:T);
          obj.infer.smooth.mu    = obj.infer.smooth.mu(:,1:T);
          obj.infer.smooth.sigma = obj.infer.smooth.sigma(1:T);
          obj.infer.smooth.G     = obj.infer.smooth.G(1:T-1);
          for pp = 1:obj.stackptr
              obj.stack{pp,1}.infer.filter.mu    = obj.stack{pp,1}.infer.filter.mu(:,1:T);
              obj.stack{pp,1}.infer.filter.sigma = obj.stack{pp,1}.infer.filter.sigma(1:T);
              obj.stack{pp,1}.infer.smooth.mu    = obj.stack{pp,1}.infer.smooth.mu(:,1:T);
              obj.stack{pp,1}.infer.smooth.sigma = obj.stack{pp,1}.infer.smooth.sigma(1:T);
              obj.stack{pp,1}.infer.smooth.G     = obj.stack{pp,1}.infer.smooth.G(1:T-1);
              obj.stack{pp,1}.yhat               = obj.stack{pp,1}.yhat(:,1:T);
          end
      end
      
      function rmNaNDimension(obj)
          curd = obj.d.y;
          newd = all(isnan(obj.y), 2);
          assert(~any(diff(newd) > 1), 'NaN Dimension in middle of object. Not supported');    % sanity check. Not difficult to change this.
          assert(obj.emiLinear, 'Removing dimensions from nonlinear emissions not supported'); % this is more difficult.
          newd = sum(~newd);
          
          if curd == newd; return; end
          
          obj.d.y  = newd;
          obj.y    = obj.y(1:newd,:);
          if ~isempty(obj.yhat); obj.yhat = obj.yhat(1:newd,:); end
          obj.par.H = obj.par.H(1:newd,:);
          obj.par.R = obj.par.R(1:newd,1:newd);
          if ~isempty(obj.par.C); obj.par.C = obj.par.C(1:newd,:); end  % control
          if ~isempty(obj.par.c); obj.par.c = obj.par.c(1:newd,:); end  % bias
          % --> DO NOT USE STACK HISTORY SINCE WILL NOT BE CONSONANT WITH
          % NEW DIMENSION... :S
      end
      
      % --- prototypes -------------------
      % inference / learning 
      [a,q]         = expLogJoint(obj, varargin); % Q(theta, theta_n) / free energy less entropy
%       obj = filterKalman(obj, bDoLLH, bDoValidation); % Kalman Filter
%       obj = filterExtended(obj, bDoLLH, bDoValidation); % Extended (EKF) Filter
%       obj = filterUnscented(obj, bDoLLH, bDoValidation, utpar); % Unscented (UKF) Filter
      D             = filter(obj, fType, bDoLLH, utpar, opts) % one of dynamics linear, other piped to NL.
      val           = tempFilterGrad(obj, fType, bDoLLH, utpar, opts, var) % DELETE ME!
      D             = getGradient(obj, par, doCheck) % get gradient of parameters
%       obj = smoothLinear(obj, bDoValidation); % RTS Smoother
%        obj = smoothExtended(obj, bDoValidation); % Extended (EKF) RTS Smoother
%       obj = smoothUnscented(obj, bDoValidation, utpar); % Unscented (UKF) RTS Smoother
      smooth(obj, fType, utpar, varargin); % one of dynamics linear, other piped to NL.
      ssid(obj, L);  % Subspace ID
      [llh, niters] = parameterLearningEM(obj, opts);
      lhHist        = parameterLearnOnline(obj, fType, opts, utpar)
      
      s             = suffStats(obj, opts);      % Sufficient statistics required for learning
      [y, covY]     = impute_y(obj, varargin);   % impute missing values into y ('filter'/'smooth', true);
      % graphical
      plotStep2D(obj, posteriorType);
      
      
      function set_y(obj, val)   % not the official setter (N/A), since accesses property d => can cause initialisation probs.
          assert(all(size(val) == [obj.d.y, obj.d.T]), 'value must be of size (%d, %d)', obj.d.y, obj.d.T);
          obj.y = val;
      end
   end
   
   methods (Access = public, Hidden=true)
       
       parameterLearningMStep(obj, updateOnly, opts); % internals for EM
       ok  = validationInference(obj, doError); % Input validation
       out = parameterHash(obj);
       
       % ---- Dynamics wrappers --------------------------
        function out = doTransition(obj, input, u)
           
            if obj.evoLinear
                  out        = obj.par.A * input;
                  if obj.hasControl(1); out = out + obj.par.B * u; end
            else
                [f,~,~,~]    = obj.functionInterfaces;
                if ~obj.hasControl(1), out = f(input);
                else, out = f(input, u); end
            end
        end
        function out = doEmission(obj, input, u, P, utpar)
            b       = zeros(obj.d.y, 1);
            if ~isempty(obj.par.c); b = obj.par.c; end
            if obj.emiLinear
                 out = obj.par.H * input + b;
                 if obj.hasControl(2); out = out + obj.par.C * u; end
            else
                [~,~,h,~]    = obj.functionInterfaces;
                if nargin == 5 && ~isempty(utpar)
                    nlpars.f = h;
                    nlpars.Q = 0;
                    if ~obj.hasControl(2); u = []; end
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

             %% ---- Save stack handlers -----------------------
      function stackPush(obj, insertion, descr)
          maxIdx = size(obj.stack, 1);
          if obj.stackptr >= maxIdx
              obj.stack = [obj.stack; cell(100,2)];
          end
          if any(strcmpi(descr, obj.getStackDescrList))
              error('Name is not unique. Please choose a different name');
          end
          obj.stackptr = obj.stackptr + 1;
          obj.stack{obj.stackptr,1} = insertion;
          obj.stack{obj.stackptr,2} = descr;
      end
      
      function stackOverwrite(obj, insertion, descr, pos)
          maxIdx = size(obj.stack, 1);
          if pos >= maxIdx
              obj.stack = [obj.stack; cell(100,2)];
          end
          obj.stackptr = max(obj.stackptr, pos);
          obj.stack{pos,1} = insertion;
          obj.stack{pos,2} = descr;
      end
      
      function stackDelete(obj, pos)
          if obj.stackptr <= 0
              fprintf('Empty save stack. Nothing to do.\n');
              return;
          end
          if nargin > 1
              idx = stackFind(obj, pos);
          else
              idx = obj.stackptr;
          end
          
          maxIdx = size(obj.stack, 1);
          obj.stack(idx,:) = [];
          obj.stackptr = obj.stackptr - 1;
            
          if obj.stackptr < maxIdx - 150
              obj.stack = obj.stack(1:(maxIdx - 100),:);
          end
      end
        
      function [contents, descr] = stackTop(obj)
          if obj.stackptr <= 0
              fprintf('Empty save stack. \n');
              contents = []; descr = []; return
          end
          contents = obj.stack{obj.stackptr,1};
          descr    = obj.stack{obj.stackptr,2};
      end
       
      function descr = getStackDescrList(obj)
          if obj.stackptr == 0
              descr = {}; return
          end
          descr = obj.stack(1:obj.stackptr,2);
      end
      
      function idx = stackFind(obj, test)
        if isnumeric(test) && isscalar(test)
            assert(test <= obj.stackptr, 'save point index requested exceeds the save stack')
            assert(test > 0, 'save point index must be positive');
            idx = test;
        elseif ischar(test)
            descr    = obj.getStackDescrList;
            idx      = find(strcmpi(test, descr));
            if isempty(idx)
                error('ds:saveNameNotOnStack','Cannot find save-point ''%s'' on save stack', test);
            end
            assert(isscalar(idx), 'Multiple matches on save stack. Should be impossible.. eep');
        else
            error('unknown search type for stack. Expected scalar-numeric, character or empty');
        end
      end
      % ------------------------------------------------
      

        
   end
   
   methods (Static)
       function val = isnummat(x)
           val = isnumeric(x) && ismatrix(x);
       end
       
       function val = isnumscal(x)
           val = isnumeric(x) && isscalar(x);
       end
       
       val = obj.parameterHash;
   end
end