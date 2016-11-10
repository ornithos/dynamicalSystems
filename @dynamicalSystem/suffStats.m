function s = suffStats(obj, opts)
% s = suffStats(obj, opts)
% Calculates the various quantities used in S�rkk�'s parameter estimation
% equations. They are typically pairwise expectations.
% opts: verbose, bIgnoreHash, bDoValidation
if nargin < 2 || isempty(opts)
    opts = struct;
end

optsDefault  = struct('verbose', true, 'bIgnoreHash', false);
opts         = utils.struct.structCoalesce(opts, optsDefault);
    
% Check for existence of Smoothed estimates
if ~opts.bIgnoreHash && ~strcmp(obj.infer.fpHash, obj.parameterHash)
    if opts.verbose; fprintf('Filter not run or parameters changed. Rerunning filter...\n'); end
    obj = obj.filter('Kalman', false, [], opts);
    obj = obj.smooth('Linear', [], opts);
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

%% calculate quantities � la S�rkk�
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

end