function prev = updateParams_bspl(obj, varargin)
    
    dy                        = obj.d.y;
    dx                        = obj.d.x;
    
    switch numel(varargin) 
        case 1
            % updateParams_bspl(obj, x)
            %   * Resize, reshape and distribute parameters (*BELOW*)
            x       = varargin{1};
            etaMask = [];
        case 2
            % updateParams_bspl(obj, x, etaMask)
            % mask out certain elements of eta
            x       = varargin{1};
            etaMask = varargin{2};
        case 3
            % updateParams_bspl(obj, eta, C, bias)
            %   * distribute parameters
            assert(nargout == 0, 'Unable to produce output for the 2-parameter specification');
            assert(all(size(varargin{1}) == size(obj.par.emiNLParams.eta)), 'eta wrong size');
            assert(all(size(varargin{2}) == [dy, dx]), 'eta wrong size');
            assert(all(size(varargin{3}) == [dy, 1]), 'bias wrong size');
            obj.par.emiNLParams.eta   = varargin{1};
            obj.par.emiNLParams.C     = varargin{2};
            obj.par.emiNLParams.bias  = varargin{3};
            
            % (----> RETURN <-----)
            return
        otherwise
            error('Do not know what to do with %d inputs. Chppse either {x} or {eta, C, bias}');
    end
    
    % update parameters
    if ~isempty(etaMask)
        szEta                 = sum(etaMask);
    else
        szEta                 = size(obj.par.emiNLParams.eta,2);
        etaMask               = true(1, szEta);
    end
    
    newEtas                   = reshape(x(1:dy*szEta), dy, szEta);
    newEtas                   = cumsum(exp(newEtas), 2);   % CHANGEME
    
    if nargout > 0
        prev = getCurrParamVector(obj, etaMask);
    end
    
    obj.par.emiNLParams.eta(:,etaMask)   = newEtas;
    
    % other params
    obj.par.emiNLParams.C     = reshape(x((dy*szEta+1):(dy*szEta+dy*dx)), dy, dx);
    obj.par.emiNLParams.bias  = reshape(x((dy*szEta+dy*dx+1):end), dy, 1);
end

function out =  getCurrParamVector(obj, etaMask)
    emiParams = obj.par.emiNLParams;
    cEta      = emiParams.eta(:, etaMask);
    out       = [cEta(:); emiParams.C(:); emiParams.bias(:)];
end