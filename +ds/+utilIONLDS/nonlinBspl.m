function [h, D] = nonlinBspl(z, pars, bspl)
    X    = pars.C * z;
    if isfield(pars.bias) && ~isempty(pars.bias)
        X = bsxfun(@plus, X, pars.bias);
    end
    
    h    = bspl.functionEval(X, pars.eta);
    
    % do derivative
    if nargout > 1
        deriv = 1;
        ltK          = X < bspl.t(1);
        gtK          = X > bspl.t(end);
        X(ltK)       = bspl.t(1);
        X(gtK)       = bspl.t(end);
        N            = size(X,1);
        D            = zeros(size(X));
        [~, m2]      = size(pars.eta);
        knots        = bspl.t;
        knots(1:(bspl.bdknots-1)) = [];
        knots((end-(bspl.bdknots-1)+1):end) = [];
        for ii = 1:N
            dB      = external.fda.bsplineM(X(ii,:), knots, bspl.k+1, deriv, 0); 
            D(ii,:) = dB*pars.eta(:,min(ii,m2));
        end
    end
end