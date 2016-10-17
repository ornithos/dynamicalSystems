classdef dynamicalSystem
   properties
      x, y, A, H, Q, R, x0, posterior=struct, d
   end
   properties (SetAccess = protected)
       inA = []
       inH = []
       inQ = []
       inR = []
       T = []
       opts = []
       hasOpts = false
       llh = [];
   end
   methods
      function obj = dynamicalSystem(x0, varargin)
         nargs = length(varargin);
         
         % find opts object (if any)
         for ii = 1:nargs
             if isa(varargin{ii}, 'dsOpts')
                 obj.opts = varargin{ii};
                 obj.hasOpts = true;
                 varargin(ii) = [];
                 nargs = nargs - 1;
             end
         end
                 
         switch nargs
             case 2
                 % input data and dimension of latent space.
                 assert(obj.isnummat(varargin{1}), 'prototype (y, d): argument one must be a numeric matrix');
                 assert(obj.isnumscal(varargin{2}), 'prototype (y, d): argument two must be a numeric scalar');
                 obj.y = varargin{1};
                 obj.d = struct('x', varargin{2}, 'y', size(obj.y, 2));
                 
                 obj   = obj.processX0(x0);
                 
                 if ~(obj.hasOpts && ~obj.opts.warnings)
                     if size(obj.y,1) > obj.d.y
                         warning('more dimensions in observations than timepoints. y is (d x T) matrix.');
                     end
                     if obj.d.y < obj.d.x
                         warning('latent space dimensionality is greater than input space!')
                     end
                 end
                 
             case {5,6}
                 % all system matrices
                 assert(obj.isnummat(varargin{1}), 'prototype (A, H, Q, R, [T/y]): argument 1 must be a numeric matrix');
                 assert(obj.isnummat(varargin{2}), 'prototype (A, H, Q, R, [T/y]): argument 2 must be a numeric matrix');
                 assert(obj.isnummat(varargin{3}), 'prototype (A, H, Q, R, [T/y]): argument 3 must be a numeric matrix');
                 assert(obj.isnummat(varargin{4}), 'prototype (A, H, Q, R, [T/y]): argument 4 must be a numeric matrix');
                 obj.inA = varargin{1};
                 obj.inH = varargin{2};
                 obj.inQ = varargin{3};
                 obj.inR = varargin{4};
                 obj.d = struct('x', size(obj.inA,1), 'y', size(obj.inH, 1));
                 assert(all(size(obj.inA)==[obj.d.x, obj.d.x]), 'matrix A is not square');
                 assert(all(size(obj.inQ)==[obj.d.x, obj.d.x]), 'matrix Q is not d.x by d.x');
                 assert(all(size(obj.inH)==[obj.d.y, obj.d.x]), 'matrix H is not d.y by d.x');
                 assert(all(size(obj.inR)==[obj.d.y, obj.d.y]), 'matrix R is not d.y by d.y');
                 obj = obj.processX0(x0);
                 
                 if obj.isnumscal(varargin{5})
                     obj.T = varargin{5};
                     obj   = obj.generateData;
                 else
                     assert(obj.isnummat(varargin{5}), 'prototype (A, H, Q, R, [T/y]): argument 5 must be a numeric scalar or matrix');
                     obj.y = varargin{5};
                     obj.T = size(obj.y, 2);
                     if ~(obj.hasOpts && ~obj.opts.warnings)
                         if size(obj.y,1) > obj.d.y
                             warning('more dimensions in observations than timepoints. y is (d x T) matrix.');
                         end
                     end
                 end
                 
             otherwise
                 error('Number of arguments (%d) does not match any prototype format!', nargs);
             
         end
         % pre populate filter/smoother for input data
         fprintf('(%s) Running smoother for input parameters...  ', datestr(now, 'HH:MM:SS'));
         obj = obj.useInputParameters;
         obj = obj.posteriorFilter;
         obj = obj.posteriorSmooth;
         obj.posterior.inFilter = obj.posterior.filter;
         obj.posterior.inSmooth = obj.posterior.smooth;
         fprintf('Complete!\n');
         obj = obj.clearParameters;
      end
      
      function obj = generateData(obj)
          obj.x       = zeros(obj.d.x, obj.T+1);
          obj.y       = zeros(obj.d.y, obj.T);
          obj.x(:,1)  = obj.x0.mu;
          
          transChol   = chol(obj.inQ);
          emissChol   = chol(obj.inR);
          for tt = 1:obj.T
              obj.x(:,tt+1) = obj.inA * obj.x(:,tt) + transChol * randn(obj.d.x,1);
              obj.y(:,tt)   = obj.inH * obj.x(:,tt+1) + emissChol * randn(obj.d.y, 1);
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
      
      % prototypes
      obj = posteriorFilter(obj, bDoLLH, bDoValidation); % Kalman Filter
      obj = posteriorSmooth(obj, bDoValidation); % RTS Smoother
      obj = validationInference(obj); % Input validation
      obj = ssid(obj, L);  % Subspace ID
      [obj, llh] = parameterLearningEM(obj, opts);
      obj = parameterLearningMStep(obj, opts);
      obj = calcLogLikelihood(obj);
      
      plotStep2D(obj, posteriorType)
   end
   
   methods (Access = private)
       function obj = processX0(obj, x0)
       % process x0 object
         if isempty(x0)
             obj.x0 = struct('mu', zeros(obj.d.x,1), 'sigma', eye(obj.d.x)*1e9);
         else
             assert(isstruct(x0) && all(ismember(fieldnames(x0), {'mu','sigma'})), ...
                 'x0 must be a struct with fields ''mu'' and ''sigma''');
             obj.x0 = x0;
         end
       end
       
       function obj = clearParameters(obj)
           obj.A = [];
           obj.H = [];
           obj.Q = [];
           obj.R = [];
       end
   end
   
   methods (Static)
       function val = isnummat(x)
           val = isnumeric(x) && ismatrix(x);
       end
       function val = isnumscal(x)
           val = isnumeric(x) && isscalar(x);
       end
       
   end
end