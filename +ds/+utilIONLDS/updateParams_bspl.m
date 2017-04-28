function prev = updateParams_bspl(obj, x, etaMask, varargin)
    
    dy                        = obj.d.y;
    dx                        = obj.d.x;
    
    logSpace                  = true;
    
    if nargin < 3
        etaMask = [];
    end
    
    if nargin >= 4
        while ~isempty(varargin)
            assert(ischar(varargin{1}), 'varargin must come in name-value pairs');
            switch varargin{1}
                case 'logSpace'
                    assert(islogical(varargin{2}) && isscalar(varargin{2}), 'logSpace must be logical scalar');
                    logSpace = varargin{2};
                otherwise
                    warning('unknown options specified: %s (I understand %s)', varargin{1}, strjoin({'logSpace'}, ','));
            end
            varargin(1:2) = [];
        end
    end
    
    
    % ________ update parameters __________________________
    
    % Eta
    if ~isempty(etaMask)
        szEta                 = sum(etaMask);
    else
        szEta                 = size(obj.par.emiNLParams.eta,2);
        etaMask               = true(1, szEta);
    end
    
    newEtas                   = reshape(x(1:dy*szEta), dy, szEta);
    
    if logSpace
        newEtas                   = exp(newEtas);
    end
    newEtas                   = cumsum(newEtas, 2);   % CHANGEME
    
    % (In case user requests output)
    if nargout > 0
        prev = ds.utilIONLDS.getCurrParamVector(obj, etaMask, logSpace);
    end
    
    % eta update
    obj.par.emiNLParams.eta(:,etaMask)   = newEtas;
    
    % other params
    obj.par.emiNLParams.C     = reshape(x((dy*szEta+1):(dy*szEta+dy*dx)), dy, dx);
    obj.par.emiNLParams.bias  = reshape(x((dy*szEta+dy*dx+1):end), dy, 1);
end

