function [f, d] = derivEmiWrapper_bspl(obj, x, varargin)
    % [f, d] = derivEmiWrapper(obj, x)
    % Wrapper for optimisation of derivative of emission parameters for *gen_sigmoid*.
    % Allows function call from vector of parameters and returns vector of
    % gradient as well as objective value.
    % Vector of parameters: x = [eta(:); C(:); bias(:)];
    
    assert(isa(obj, 'ds.ionlds'),  'input object is not a valid IONLDS object');
    dy    = obj.d.y;
    dx    = obj.d.x;
    nKnot = 6;
    assert(isnumeric(x) && numel(x) == dy*nKnot + dy*(dx+1), 'parameter vector is not the correct size');
    
    % save current parameters
    eta   = obj.par.emiNLParams.eta;
    C     = obj.par.emiNLParams.C;
    bias  = obj.par.emiNLParams.bias;
    
    % update parameters from input
    ds.utilIONLDS.updateParams_bspl(obj, x);
    
    % get function val and gradient
    if nargout > 1
        [f, D] = obj.expLogJoint_bspl(varargin{:});
        
        % reorganise output
        d                = zeros(size(x));
        consec           = 1:obj.d.y;

        d(dy*0 + consec) = D.m;
        d(dy*1 + consec) = D.M;
        d(dy*2 + consec) = D.nu;
        d(dy*3 + consec) = D.gamma;
    %     d(dy*4+1:end) = 0;
        d((dy*4 + 1):end)= D.C(:);
        
        d                = -d;  % see below (negation)
    else
        f = obj.expLogJoint_bspl(varargin{:});
    end
    
    % reset object parameters (is a handle object)
    obj.par.emiNLParams.eta   = eta;
    obj.par.emiNLParams.C     = C;
    obj.par.emiNLParams.bias  = bias;
    
    % we are trying to maximise, but generic optimisation functions
    % minimise, so we negate the function value
    f                = -f;
end