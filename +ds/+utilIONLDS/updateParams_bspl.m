function updateParams_bspl(obj, varargin)
    
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
        szEta                     = sum(etaMask);
        obj.par.emiNLParams.eta(:,etaMask)   = reshape(x(1:dy*szEta), dy, szEta);
    else
        eta                       = obj.par.emiNLParams.eta;
        szEta                     = size(eta,2);
        obj.par.emiNLParams.eta   = reshape(x(1:dy*szEta), dy, szEta);
    end
    obj.par.emiNLParams.C     = reshape(x((dy*szEta+1):(dy*szEta+dy*dx)), dy, dx);
    obj.par.emiNLParams.bias  = reshape(x((dy*szEta+dy*dx+1):end), dy, 1);
end