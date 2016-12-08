classdef ionlds < ds.dynamicalSystem
   % ionlds: create input-output non-linear dynamical system object à la
   % (Georgatzis et al. 2016) inheriting from the general dynamicalSystem
   % class. This differs only in the EM procedure to learn the nonlinear
   % function parameters. At present dynamicalSystem cannot learn these,
   % and it is difficult to write a performant general purpose optimiser
   % for these (although eg. autodiff could be done).
   %
   % ionlds(dimx, dimy, 'evolution', {A || f || []}, {Df, Q}
   %                 'emission', {H || h || []}, {Dh, R},
   %                 'data', {y || T}, 'x0' {x0mu, x0cov}, 
   %                 ('xtrue', x),
   %                 ('control', u, (C), (D)),
   %                 opts)
   %
   % Please refer to dynamicalSystem for documentation on the constructor.
   % ionlds uses the same syntax but enforces that the latent dynamics are
   % linear.
   
   methods
      function obj = ionlds(varargin)
         % CONSTRUCTOR
         obj@ds.dynamicalSystem(varargin{:});  % do superclass constructor
         assert(obj.evoLinear, 'transition model must be linear in IONLDS');
         assert(~obj.emiLinear, 'emissions must be non-linear in IONLDS');
         assert(obj.hasControl(1), 'transition model must include inputs in IONLDS');
         assert(~obj.hasControl(2), 'emission model must not include inputs in IONLDS');
      end
      
      % modified superclass methods
      [a,q]              = expLogJoint(obj); % Q(theta, theta_n) / free energy less entropy
      D                  = getGradient(obj, par, doCheck) % get gradient of parameters
      [llh, niters]      = parameterLearningEM(obj, opts); % do learning
   end
   
   methods (Access = protected, Hidden=true)
       parameterLearningMStep(obj, updateOnly, opts); % internals for EM
   end
end