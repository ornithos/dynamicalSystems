function [f, d] = derivEmiWrapper_bspl(obj, x, varargin)
    % [f, d] = derivEmiWrapper(obj, x)
    % Wrapper for optimisation of derivative of emission parameters for *gen_sigmoid*.
    % Allows function call from vector of parameters and returns vector of
    % gradient as well as objective value.
    % Vector of parameters: x = [eta(:); C(:); bias(:)];
    
    assert(isa(obj, 'ds.ionlds'),  'input object is not a valid IONLDS object');
    dy    = obj.d.y;
    dx    = obj.d.x;
    
    nvgin = numel(varargin);
    
    % process varargin manually for speed
    opts.fixBias2 = false;
    opts.utpar    = struct;
    opts.etaMask  = [];
    if nvgin > 3 && ischar(varargin{1}) 
        while ~isempty(varargin)
            assert(ischar(varargin{1}), 'varargin must come in name-value pairs');
            switch varargin{1}
                case 'fixBias2'
                    assert(islogical(varargin{2}) && isscalar(varargin{2}), 'fixBias2 must be logical scalar');
                    opts.fixBias2 = varargin{2};
                case 'utpar'
                    assert(isstruct(varargin{2}), 'utpars must be a struct');
                    if all(ismember({'alpha', 'beta', 'kappa'}, fieldnames(varargin{2})))
                        opts.utpars = varargin{2};
                    elseif ~isempty(fieldnames(varargin{2}))
                        warning('Need alpha, beta, kappa fieldnames in utpars. Not found so ignoring');
                    end
                case 'etaMask'
                    assert(numel(varargin{2}) == size(obj.par.emiNLParams.eta,2), 'etaMask is not conformable to eta');
                    opts.etaMask = varargin{2};
                otherwise
                    warning('unknown options specified: %s (I understand %s)', varargin{1}, strjoin(opts.fieldnames, ','));
            end
            varargin(1:2) = [];
        end
    end
    
    nKnot = sum(opts.etaMask);
    assert(isnumeric(x) && numel(x) == dy*nKnot + dy*(dx+1), 'parameter vector is not the correct size');
    
    fixBias2    = opts.fixBias2;
    utpar       = opts.utpar;
    opts        = rmfield(opts, 'fixBias2');  % for use in bsplineGrad
    opts        = rmfield(opts, 'utpar');
    %% Function body
    
    % save current parameters
    eta   = obj.par.emiNLParams.eta;
    C     = obj.par.emiNLParams.C;
    bias  = obj.par.emiNLParams.bias;
    
    % update parameters from input
    ds.utilIONLDS.updateParams_bspl(obj, x, opts.etaMask);
    
    % replace bias so as not to mess with finite difference (if checking)
    if fixBias2
        obj.par.emiNLParams.bias    = bias;
    end
    
    % get function val and gradient
    if nargout > 1
        [f, D] = ds.utilIONLDS.bsplineGrad(obj, utpar, opts);
        if ~fixBias2
            d                = [D.eta(:); D.C(:); D.bias(:)];
        else
            d                = [D.eta(:); D.C(:); zeros(size(D.bias(:)))];
        end
        % maximisation ----> minimisation
        f                = -f;
        d                = -d;  % see below (negation)
    else
        f                = ds.utilIONLDS.bsplineGrad(obj, utpar, 'gradient', false);
        
        % maximisation ----> minimisation
        f                = -f;
%         [~,M2] = ds.utilIONLDS.utTransform_ymHx_bspl(obj);
%         f = -0.5*trace(inv(obj.par.R)*M2);
    end
    
    % reset object parameters (is a handle object)
    ds.utilIONLDS.updateParams_bspl(obj, eta, C, bias);
    
end