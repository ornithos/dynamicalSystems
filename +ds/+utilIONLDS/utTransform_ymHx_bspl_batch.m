function [ymHx, outerprod] = utTransform_ymHx_bspl_batch(obj, alpha, beta, kappa)
    % wrapper for utTransform_ymHx_bspl which accepts batch object inputs
    assert(isa(obj, 'ds.dynamicalSystemBatch'), 'non-batch input to batch function!');
    
    if nargin == 1
        alpha = 1;
        beta  = 0;
        kappa = 0;
    end
        
    dsObjects = cell(obj.d.n,1);
    for nn = 1:obj.d.n; dsObjects{nn} = obj.extractSingle(nn, false); end
    N         = obj.d.n;
    ymHx      = zeros(obj.d.y, N);
    outerprod = zeros(obj.d.y, obj.d.y);
    
    if nargout == 1    
        for nn = 1:N
            ymHx(:,nn) = ds.utilIONLDS.utTransform_ymHx_bspl(dsObjects{nn}, alpha, beta, kappa);
        end
    else
        for nn = 1:N
            [yy,rr]    = ds.utilIONLDS.utTransform_ymHx_bspl(dsObjects{nn}, alpha, beta, kappa);
            ymHx(:,nn) = mean(yy,2);
            outerprod  = outerprod + rr;
        end
    end
    
    ymHx = sum(ymHx,2)./(sum(ymHx~=0,2));
end