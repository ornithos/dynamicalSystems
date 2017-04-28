classdef ionldsBatch < ds.dynamicalSystemBatch & ds.ionlds
   % ionldsBatch: create input-output non-linear dynamical system object à
   % la (Georgatzis et al. 2016) inheriting from the general dynamicalSystem
   % class. This is the batch learning variant:
   %           ionldsBatch <-> ionlds 
   %            in the same way that
   %  dynamicalSystemBatch <-> dynamicalSystem.
   %
   % ionldsBatch(dimx, dimy, 'evolution', {A || f || []}, {Df, Q}
   %                 'emission', {H || h || []}, {Dh, R},
   %                 'data', {y || T}, 'x0' {x0mu, x0cov}, 
   %                 ('xtrue', x),
   %                 ('control', u, (C), (D)),
   %                 opts)
   %
   % Please refer to dynamicalSystemBatch for documentation on the constructor.
   % ionldsBatch uses the same syntax but enforces that the latent dynamics are
   % linear.
   
   methods
      function obj = ionldsBatch(varargin)
         % CONSTRUCTOR
         obj@ds.dynamicalSystemBatch(varargin{:});  % do superclass constructor 
         obj@ds.ionlds('nullConstruct');
         
         assert(obj.evoLinear, 'transition model must be linear in IONLDS');
         assert(~obj.emiLinear, 'emissions must be non-linear in IONLDS');
         assert(obj.hasControl(1), 'transition model must include inputs in IONLDS');
         assert(~obj.hasControl(2), 'emission model must not include inputs in IONLDS');
      end
      
      % Make a copy of a handle object.
      function new = copy(this)
          % Instantiate new object of the same class.
          curWarns           = this.opts.warnings;
          this.opts.warnings = false;
          new                = ds.ionldsBatch(this);
          this.opts.warnings = curWarns;
      end
      
      % cannot initialise with ionlds superclass, so hacking the ionlds
      % methods onto ionldsBatch:
%       function [llh, niters]      = parameterLearningEM_bspl(obj, opts)
%           [llh, niters] = parameterLearningEM_bspl@ds.ionlds(obj, opts);
%       end
%       function [llh, niters]      = expLogJoint_bspl(obj, opts)
%           [llh, niters] = ds.ionlds.expLogJoint_bspl(obj, opts);
%       end

   end
      
end