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
    
    % process varargin manually for speed (? not sure this is still worth
    % it -- initially were only a couple of options)
    fixBias2 = false;
    utpar    = struct;
    etaUB    = Inf(dy,1);
    opts.etaMask  = [];
    
    if nvgin > 3 && ischar(varargin{1}) 
        while ~isempty(varargin)
            assert(ischar(varargin{1}), 'varargin must come in name-value pairs');
            switch varargin{1}
                case 'fixBias2'
                    assert(islogical(varargin{2}) && isscalar(varargin{2}), 'fixBias2 must be logical scalar');
                    fixBias2 = varargin{2};
                case 'utpar'
                    assert(isstruct(varargin{2}), 'utpars must be a struct');
                    if all(ismember({'alpha', 'beta', 'kappa'}, fieldnames(varargin{2})))
                        utpar = varargin{2};
                    elseif ~isempty(fieldnames(varargin{2}))
                        warning('Need alpha, beta, kappa fieldnames in utpars. Not found so ignoring');
                    end
                case 'etaMask'
                    assert(numel(varargin{2}) == size(obj.par.emiNLParams.eta,2), 'etaMask is not conformable to eta');
                    opts.etaMask = varargin{2};
                case 'etaUB'
                    if ~isempty(varargin{2})
                        assert(numel(varargin{2}) == size(obj.par.emiNLParams.eta,1), 'etaUB must have d dimensions');
                        etaUB = varargin{2};
                    end
                case 'bfgsSpline'
                    if ~isempty(varargin{2})
                        assert(isscalar(varargin{2}) && islogical(varargin{2}), 'bfgsSpline must be scalar logical');
                        bfgsSpline = varargin{2};
                    end
                otherwise
                    warning('unknown options specified: %s (I understand %s)', varargin{1}, strjoin(opts.fieldnames, ','));
            end
            varargin(1:2) = [];
        end
    end
    
    nKnot = sum(opts.etaMask);
    assert(isnumeric(x) && numel(x) == dy*nKnot + dy*(dx+1), 'parameter vector is not the correct size');

    %% Function body
    
    % save current parameters
%     eta   = obj.par.emiNLParams.eta;
%     C     = obj.par.emiNLParams.C;
%     bias  = obj.par.emiNLParams.bias;
%     
%     % update parameters from input
%     ds.utilIONLDS.updateParams_bspl(obj, x, opts.etaMask);
    
    % replace bias so as not to mess with finite difference (if checking)
    if fixBias2
        obj.par.emiNLParams.bias    = bias;
    end
    
    % get function val and gradient
    if nargout > 1
        [f, D] = ds.utilIONLDS.bsplineGradMono(obj, x, utpar, opts);
        if ~fixBias2
            d                = [D.eta(:); D.C(:); D.bias(:)];
        else
            d                = [D.eta(:); D.C(:); zeros(size(D.bias(:)))];
        end
        
        % (ie do not optimise spline coefficients)
        if ~bfgsSpline
            d(1:nKnot*dy)    = 0;
        end
%         
%         violateUB        = bsxfun(@gt, obj.par.emiNLParams.eta(:,opts.etaMask), etaUB);   % <-- this is nonsense since obj pars changed *INSIDE* grad fn now...
%         violateUB        = [violateUB(:); false(dy*(dx+1),1)] & d > 0;
%         if any(any(violateUB))
%            d(violateUB)     = 0;
%         end
        
        % maximisation ----> minimisation
        f                = -f;
        d                = -d;  % see below (negation)
    else
        f                = ds.utilIONLDS.bsplineGradMono(obj, x, utpar, 'gradient', false);
        
        % maximisation ----> minimisation
        f                = -f;
%         [~,M2] = ds.utilIONLDS.utTransform_ymHx_bspl(obj);
%         f = -0.5*trace(inv(obj.par.R)*M2);
    end
    
    % reset object parameters (is a handle object)
%     ds.utilIONLDS.updateParams_bspl(obj, eta, C, bias);
    
    
end