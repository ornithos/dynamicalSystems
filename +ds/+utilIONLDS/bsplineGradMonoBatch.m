function [f, grad, more] = bsplineGradMonoBatch(obj, x, utpar, varargin)
    % wrapper for bsplineGradMono which accepts batch object inputs
    assert(isa(obj, 'ds.dynamicalSystemBatch'), 'non-batch input to batch function!');
    
    dsObjects = cell(obj.d.n,1);
    for nn = 1:obj.d.n; dsObjects{nn} = obj.extractSingle(nn, false); end
    N = obj.d.n;
    
    if nargout == 1
        f = 0;
        for nn = 1:N
            f = f + ds.utilIONLDS.bsplineGradMono(dsObjects{nn}, x, utpar, varargin{:});
        end
    end
    if nargout > 1
        f = 0;
        g = struct('eta', 0, 'C', 0, 'bias', 0);
        more = struct('v', 0, 'K', 0, 'Hx', 'N/A for batch => to save memory');
        for nn = 1:N
            [ff, gg, mm] = ds.utilIONLDS.bsplineGradMono(dsObjects{nn}, x, utpar, varargin{:});
            f = f + ff;
            g.eta   = g.eta + gg.eta;
            g.C     = g.C   + gg.C;
            g.bias  = g.bias+ gg.bias;
            more.v  = more.v + mm.v;
            more.K  = more.K + mm.K;
%             more.Hx = more.Hx + mm.Hx;
        end
        grad = g;
    end
end
        
