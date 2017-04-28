function [m, sigma, covxy] = unscentedTransform(fn, mu, sigma, alpha, beta, kappa)
    % unscentedTransform(fn, mu, sigma)
    % Unscented transform as per (e.g.) Julier & Uhlmann 1996. Propagates
    % uncertainty described by a Guassian distribution with mean mu and
    % covariance sigma through a function and calculates closest
    % transformed distribution in the Gaussian family. This is done using
    % deterministic sigma point methods similar to Quadrature methods.
    %
    % Since many calls will be made to this function, no error checking is
    % present: this is presumed to be done before the function.
    %
    % INPUTS:
    % fn     - a function handle capable of ingesting a matrix (where each
    %          datapoint corresponds to a column) and outputting a matrix
    %          (where each column corresponds to the mapped datapoint).
    % mu     - The mean of the initial Gaussian distribution. Must be of
    %          same dimension as number of rows in the matrix described
    %          above.
    % sigma  - The covariance matrix of the Gaussian distribution.
    %
    % OPTIONAL:
    % All the below parameters may be specified or none of them.
    % alpha  - scaling parameter (related to how clustered the sigma points
    %          are around the mean.) Suggest 10^-3
    % beta   - adjustment made for kurtosis of transformed distn. Suggest 2
    % kappa  - original scaling parameter of UT. Suggest (3 - dim(mu))
    %
    % OUTPUTS:
    % m      -  Mean of the closest transformed distribution according to
    %          the Unscented transform.
    % sigma  -  Covariance matrix of the same distribution.
    %
    
    n        = size(mu, 1);
    if nargin == 3
        alpha = 1;
        beta  = 0;
        kappa = 3-n;
    end
    lambda   = alpha^2 * (n + kappa) - n;
    scl      = sqrt(n + lambda);
    
    % covariance spread
    S        = scl.*chol(sigma, 'lower');
    
    % sigma points
    spts     = [zeros(n,1), S, -S];
    spts     = bsxfun(@plus, spts, mu);
    ysp      = fn(spts);
    
    % weights
    Wm       = ones(2*n+1, 1)./(2*(n + lambda));
    Wm(1)    = lambda/(n+lambda);
    Wc       = Wm;
    Wc(1)    = Wc(1) + (1 - alpha^2 + beta);
%     sWc      = sqrt(Wc);
    
    % calculate transformed
    m        = ysp * Wm;
    delta    = bsxfun(@minus, ysp, m);
    sigma    = bsxfun(@times, delta, Wc') * delta';
    
    if nargout  == 3
        cen_spts = bsxfun(@minus, spts, mu);
        covxy    = bsxfun(@times, cen_spts, Wc') * delta';
    end
    
end