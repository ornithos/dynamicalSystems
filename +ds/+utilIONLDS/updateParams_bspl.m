function updateParams_bspl(obj, x)
    dy                        = obj.d.y;
    dx                        = obj.d.x;
    eta                       = obj.par.emiNLParams.eta;
    szEta                     = size(eta,2) -4;    % don't optimise (*2 into) endpoints
    
    % update parameters
    obj.par.emiNLParams.eta   = [eta(:, 1:2), reshape(x(1:dy*szEta), dy, szEta), eta(:, end-1:end)];
    obj.par.emiNLParams.C     = reshape(x((dy*szEta+1):(dy*szEta+dy*dx)), dy, dx);
    obj.par.emiNLParams.bias  = reshape(x((dy*szEta+dy*dx+1):end), dy, 1);
end