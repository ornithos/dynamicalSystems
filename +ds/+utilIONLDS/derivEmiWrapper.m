function [f, d] = derivEmiWrapper(obj, x)
    % [f, d] = derivEmiWrapper(obj, x)
    % Wrapper for optimisation of derivative of emission parameters. Allows
    % function call from vector of parameters and returns vector of
    % gradient as well as objective value.
    % Vector of parameters: x = [eta(:); C(:)];
    
    assert(isa(obj, 'ds.ionlds'),  'input object is not a valid IONLDS object');
    dy    = obj.d.y;
    dx    = obj.d.x;
    assert(isnumeric(x) && numel(x) == dy*4 + dy*dx, 'parameter vector is not the correct size');
    
    % save current parameters
    eta   = obj.par.emiNLParams.eta;
    C     = obj.par.emiNLParams.C;
    szEta = size(eta,2);
    
    % update parameters
    obj.par.emiNLParams.eta   = reshape(x(1:dy*szEta), dy, szEta);
    obj.par.emiNLParams.C     = reshape(x((dy*szEta+1):end), dy, dx);
    
    % get function val and gradient
    if nargout > 1
        [f, D] = obj.expLogJoint;
        
        % reorganise output
        d                = zeros(size(x));
        consec           = 1:obj.d.y;

        d(dy*0 + consec) = D.m;
        d(dy*1 + consec) = D.M;
        d(dy*2 + consec) = D.gamma;
        d(dy*3 + consec) = D.b;
    %     d(dy*4+1:end) = 0;
        d((dy*4 + 1):end)= D.C(:);
        
        d                = -d;  % see below (negation)
    else
        f = obj.expLogJoint;
    end
    
    % reset object parameters (is a handle object)
    obj.par.emiNLParams.eta   = eta;
    obj.par.emiNLParams.C     = C;
    
    % we are trying to maximise, but generic optimisation functions
    % minimise, so we negate the function value
    f                = -f;
end