classdef dynamicalSystem
   % dynamicalSystem: create linear/non-linear dynamical system object with
   % access to filtering, smoothing and learning methods, and visualisations
   % for comparison of methods.
   %
   % dynamicalSystem(dimx, dimy, 'evolution', {A || f || []}, {Df, Q}
   %                 'emission', {H || h || []}, {Dh, R},
   %                 'data', {y || T}, 'x0' {x0mu, x0cov}, opts)
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
   %
   % *** EMISSION ***
   % - the first argument is the same as that of the evolution parameters,
   % corresponding to the emission linear/nonlinear function.
   % - a couple of optional arguments. Dh = the derivative (Hessian) of h,
   % the emission noise covariance R (can be learned if not specified).
   %
   % *** DATA ***
   % - either a data matrix y corresponding to the observed data (each
   % column is an observation) at time 1, ..., t; or the number of
   % observations that should be generated from the given dynamics. The
   % latter may only be specified if all transitions and covariances have
   % also been specified.
   %
   % % ***** x0 *****
   % - x0mu: the mean of the prior over the latent space (if not given
   % assumed 0)
   % - x0cov: the covariance of the prior over the latent space (required).
   %
   % *** OPTS ***
   % - A struct containing fields: warnings, verbose
   
   properties
      x, y, f, Df, h, Dh, x0, stack = cell(100,2), d
   end
   properties (SetAccess = protected)
       opts = []
       evoLinear = []
       evoNLhasParams = false
       emiLinear = []
       emiNLhasParams = false
       stackptr = 0
       % ----------------- Changeable:
       A = []
       H = []
       Q = []
       R = []
       llh = []
       evoNLParams = struct
       emiNLParams = struct
       filter = []
       fpHash = []
       smooth = []
   end
   methods
      function obj = dynamicalSystem(varargin)
         % CONSTRUCTOR
         obj = obj.processInputArgs(varargin);
         
         % pre populate filter/smoother for input data
         if obj.validationInference(false)
             fprintf('(%s) Running smoother for input parameters...  ', datestr(now, 'HH:MM:SS'));
             obj = obj.filterKalman;
             obj = obj.smoothLinear;
             fprintf('Complete!\n');
             obj = obj.save('initialised-run');
         else
             if obj.opts.warnings
                 warning('dynamicalSystem cannot perform initial inference since not all parameters specified.');
             end
         end
         
      end
      
      function obj = generateData(obj)
          obj.x       = zeros(obj.d.x, obj.d.T+1);
          obj.y       = zeros(obj.d.y, obj.d.T);
          obj.x(:,1)  = obj.x0.mu;
          
          transChol   = chol(obj.Q);
          emissChol   = chol(obj.R);
          for tt = 1:obj.d.T
              obj.x(:,tt+1) = obj.doTransition(obj.x(:,tt)) + transChol * randn(obj.d.x,1);
              obj.y(:,tt)   = obj.doEmission(obj.x(:,tt+1)) + emissChol * randn(obj.d.y,1);
          end
          obj.x       = obj.x(:,2:end);
      end
        
      function obj = useInputParameters(obj)
          assert(~isempty(obj.inA), 'Cannot use input parameters: no parameters given');
          obj.A = obj.inA;
          obj.H = obj.inH;
          obj.Q = obj.inQ;
          obj.R = obj.inR;
      end
      
      
      % ---- Save stack handlers -----------------------
      function obj = stackPush(obj, insertion, descr)
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
        
      function obj = stackDelete(obj)
          if obj.stackptr <= 0
              fprintf('Empty save stack. Nothing to do.\n');
              return;
          end
          maxIdx = size(obj.stack, 1);
          obj.stack(obj.stackptr,:) = {[],[]};
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
      % ------------------------------------------------
        
      % ----- save / stack aliases ------------------------
      function obj = save(obj, descr)
          assert(nargin == 2, 'please provide a description');
          if obj.fpHash ~= obj.parameterHash
              if obj.opts.warnings; warning('posterior does not match current parameters'); end
              if obj.evoLinear && obj.emiLinear
                  fprintf('(%s) Running posterior...\n', datestr(now, 'HH:MM:SS'));
                  obj = obj.filterKalman;
                  obj = obj.smoothLinear;
              else
                  error('Unable to save: posterior does not match parameters. Please run a posterior algm');
              end
          end
          
          insertion = struct;
          insertion.A = obj.A;
          insertion.H = obj.H;
          insertion.Q = obj.Q;
          insertion.R = obj.R;
          insertion.llh = obj.llh;
          insertion.evoNLParams = obj.evoNLParams;
          insertion.emiNLParams = obj.emiNLParams;
          insertion.filter = obj.filter;
          insertion.smooth = obj.smooth;
          insertion.fpHash = obj.fpHash;
          obj = stackPush(obj, insertion, descr);
      end
      function descr = savedList(obj)
          descr = strjoin(obj.getStackDescrList, ',\n');
      end
      
      % --- prototypes -------------------
      % inference / learning 
      obj = filterKalman(obj, bDoLLH, bDoValidation); % Kalman Filter
      obj = smoothLinear(obj, bDoValidation); % RTS Smoother
      obj = validationInference(obj, doError); % Input validation
      obj = ssid(obj, L);  % Subspace ID
      [obj, llh] = parameterLearningEM(obj, opts);
      obj = parameterLearningMStep(obj, verbose, updateOnly);
      obj = calcLogLikelihood(obj);
      
      % graphical
      plotStep2D(obj, posteriorType)
   end
   
   methods (Access = private)

       % ---- Dynamics wrappers --------------------------
        function out = doTransition(obj, input)
            if obj.evoLinear
                  out = obj.A * input;
            else
                if ~obj.evoNLhasParams
                    out = obj.f(input);
                else
                    out = obj.f(input, obj.evoNLParams);
                end
            end
        end
        function out = doEmission(obj, input)
            if obj.emiLinear
                  out = obj.H * input;
            else
                if ~obj.emiNLhasParams
                    out = obj.f(input);
                else
                    out = obj.f(input, obj.emiNLParams);
                end
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