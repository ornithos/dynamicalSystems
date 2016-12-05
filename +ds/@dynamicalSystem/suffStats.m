function s = suffStats(obj, opts)
% s = suffStats(obj, opts)
% Calculates the various quantities used in Särkkä's parameter estimation
% equations. They are typically pairwise expectations.
% opts: verbose, bIgnoreHash, bDoValidation
if nargin < 2 || isempty(opts)
    opts = struct;
end

optsDefault  = struct('verbose', true, 'bIgnoreHash', false);
opts         = utils.base.parse_argumentlist(opts, optsDefault);
    
% Check for existence of Smoothed estimates
if ~opts.bIgnoreHash && obj.parametersChanged
    if opts.verbose; fprintf('Filter not run for current params. Rerunning filter/smoother...\n'); end
    obj.filter('Kalman', false, [], opts);
    obj.smooth('Linear', [], opts);
end
if isempty(obj.infer.sType)
    if opts.verbose; fprintf('Smoother not run for current params. Running smoother...\n'); end
    obj.smooth('Linear', [], opts);
end

% concatenate with x0 to provide estimates of x from 0 to T
mu    = [obj.infer.smooth.x0.mu, obj.infer.smooth.mu];
sigma = vertcat(obj.infer.smooth.x0.sigma, obj.infer.smooth.sigma);
G     = vertcat(obj.infer.smooth.x0.G, obj.infer.smooth.G);
y     = obj.y;
T     = obj.d.T;

% mu    indexes from 0:T     (T+1 elements)
% sigma indexes from 0:T     (T+1 elements)
% G     indexes from 0:(T-1) (T   elements)

%% calculate quantities à la Särkkä
almostAllP = sum(cat(3,sigma{2:T}),3);
almostAllM = mu(:,2:T) * mu(:,2:T)';
s       = struct;
s.SIGMA = (almostAllP + sigma{T+1} + almostAllM + mu(:,T+1) * mu(:,T+1)')./T;
s.PHI   = (almostAllP + sigma{1} + almostAllM + mu(:,1) * mu(:,1)')./T;
s.B     = (y * mu(:,2:T+1)')./T;
s.D     = (y * y')./T;

s.C     = mu(:,2:T+1) * mu(:,1:T)';
for tt = 1:T
    s.C = s.C + sigma{tt+1} * G{tt}';
end
s.C     = s.C./T;

s.infer = struct('mu', mu); s.infer.P = sigma; s.infer.G = G;   % hack to stop MATLAB making a struct array
end