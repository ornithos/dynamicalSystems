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
       par = struct('x0',struct,'A',[],'H',[],'Q',[],'R',[], ...
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
         
         if numel(varargin) == 1 && isa(varargin{1}, 'ds.dynamicalSystem')
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
      
      function generateData(obj)
          obj.x       = zeros(obj.d.x, obj.d.T+1);
          obj.y       = zeros(obj.d.y, obj.d.T);
          obj.x(:,1)  = obj.par.x0.mu;
          
          transChol   = chol(obj.par.Q);
          emissChol   = chol(obj.par.R);
          for tt = 1:obj.d.T
              u_t = [];
              if any(obj.hasControl); u_t = obj.u(:,tt); end
              obj.x(:,tt+1) = obj.doTransition(obj.x(:,tt), u_t) + transChol * randn(obj.d.x,1);
              obj.yhat(:,tt)= obj.doEmission(obj.x(:,tt+1), u_t);
              obj.y(:,tt)   = obj.yhat(:,tt) + emissChol * randn(obj.d.y,1);
              
              %plot(obj.par.H(1,:) * obj.x(:,1:tt-1), obj.par.H(2,:) * obj.x(:,1:tt-1)); hold on; plot(obj.y(1,1:tt-1), obj.y(2,1:tt-1)); hold off;
              %pause
          end
          obj.x       = obj.x(:,2:end);
      end
      
      % Make a copy of a handle object.
      function new = copy(this)
          % Instantiate new object of the same class.
          curWarns           = this.opts.warnings;
          this.opts.warnings = false;
          new                = ds.dynamicalSystem(this);
          this.opts.warnings = curWarns;
      end
      
      function getFittedValues(obj)
          obj.ensureInference('FITVALS', 'smooth');
          fitted = zeros(obj.d.y, obj.d.T);
          for tt = 1:obj.d.T
              u_t = [];
              if any(obj.hasControl); u_t = obj.u(:,tt); end
              fitted(:,tt) = obj.doEmission(obj.infer.smooth.mu(:,tt), u_t);
          end
          obj.yhat = fitted;
      end
      
      function fitted = getPredictedValues(obj, nlookahead)
          if nargin < 2; nlookahead = 0; end
          obj.ensureInference('PREDVALS', 'filter');
          fitted = zeros(obj.d.y, obj.d.T - nlookahead);
          for tt = 1:obj.d.T - nlookahead
              x_t = obj.infer.filter.mu(:,tt);
              u_t = [];
              if any(obj.hasControl); u_t = obj.u(:,tt); end
              for jj = 1:nlookahead
                  x_t = obj.doTransition(x_t, u_t);
                  u_t = [];
                  if any(obj.hasControl); u_t = obj.u(:,tt+jj); end
              end
              fitted(:,tt) = obj.doEmission(x_t, u_t);
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
      function [fe, Dfe, he, Dhe] = functionInterfaces(obj)
        % non-linear wrappers for common interface to fns
        if obj.evoNLhasParams
            if obj.hasControl(1)
                fe = @(x,u) obj.par.f(x, u, obj.par.evoNLParams);
                Dfe = @(x,u) obj.par.Df(x, u, obj.par.evoNLParams);
            else
                fe = @(x,u) obj.par.f(x, obj.par.evoNLParams);
                Dfe = @(x,u) obj.par.Df(x, obj.par.evoNLParams);
            end
        else
            if obj.hasControl(1)
                fe = @(x,u) obj.par.f(x,u);
                Dfe = @(x,u) obj.par.Df(x,u);
            else
                fe = @(x,u) obj.par.f(x);
                Dfe = @(x,u) obj.par.Df(x);
            end
        end
        if obj.emiNLhasParams
            if obj.hasControl(2)
                he = @(x,u) obj.par.h(x, u, obj.par.emiNLParams);
                Dhe = @(x,u) obj.par.Dh(x, u, obj.par.emiNLParams);
            else
                he = @(x) obj.par.h(x, obj.par.emiNLParams);
                Dhe = @(x) obj.par.Dh(x, obj.par.emiNLParams);
            end
        else
            if obj.hasControl(2)
                he = @(x,u) obj.par.h(x,u);
                Dhe = @(x,u) obj.par.Dh(x,u);
            else
                he = @(x,u) obj.par.h(x);
                Dhe = @(x,u) obj.par.Dh(x);
            end
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
          if ~strcmp(obj.infer.fpHash, obj.parameterHash)
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
          if nargin < 2 || isempty(fType); fType = []; end
          tmpobj = obj.copy;
          tmpobj.filter(fType, true, utpar, opts);
          llh    = tmpobj.infer.llh;
      end
      
      % --- prototypes -------------------
      % inference / learning 
      [a,q]         = expLogJoint(obj, varargin); % Q(theta, theta_n) / free energy less entropy
%       obj = filterKalman(obj, bDoLLH, bDoValidation); % Kalman Filter
%       obj = filterExtended(obj, bDoLLH, bDoValidation); % Extended (EKF) Filter
%       obj = filterUnscented(obj, bDoLLH, bDoValidation, utpar); % Unscented (UKF) Filter
      filter(obj, fType, bDoLLH, utpar, opts) % one of dynamics linear, other piped to NL.
      D             = getGradient(obj, par, doCheck) % get gradient of parameters
%       obj = smoothLinear(obj, bDoValidation); % RTS Smoother
%        obj = smoothExtended(obj, bDoValidation); % Extended (EKF) RTS Smoother
%       obj = smoothUnscented(obj, bDoValidation, utpar); % Unscented (UKF) RTS Smoother
      smooth(obj, fType, utpar, opts) % one of dynamics linear, other piped to NL.
      ssid(obj, L);  % Subspace ID
      [llh, niters] = parameterLearningEM(obj, opts);
      lhHist        = parameterLearnOnline(obj, fType, opts, utpar)
      
      s             = suffStats(obj, opts);  % Sufficient statistics required for learning
      % graphical
      plotStep2D(obj, posteriorType)
   end
   
   methods (Access = protected, Hidden=true)
       parameterLearningMStep(obj, updateOnly, opts); % internals for EM
       ok  = validationInference(obj, doError); % Input validation
       out = parameterHash(obj);
   end
   
   methods (Access = private)

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
        function out = doEmission(obj, input, u)
            if obj.emiLinear
                 out = obj.par.H * input;
                 if obj.hasControl(2); out = out + obj.par.C * u; end
            else
                [~,~,h,~]    = obj.functionInterfaces;
                if ~obj.hasControl(2), out = h(input);
                else, out = h(input, u); end
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