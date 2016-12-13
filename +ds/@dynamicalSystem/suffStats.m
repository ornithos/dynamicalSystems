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

if any(obj.hasControl)
    U = [zeros(obj.d.u, 1), obj.u];
end

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

% ----- Control quantities -----------------------------------
if obj.hasControl(1)
    s.X_Um1   = (mu(:,2:T+1) * U(:,1:T)')./T;
    s.Xm1_U   = (mu(:,1:T)   * U(:,2:T+1)')./T;
    s.U_Um1   = (U(:,2:T+1)  * U(:,1:T)')./T;
end

if obj.hasControl(2)
    s.YU      = (y * U(:,2:T+1)')./T;
end

if any(obj.hasControl)
    s.Um1Um1  = (U(:,1:T) * U(:,1:T)')./T;
    s.Xm1Um1  = (mu(:,1:T) * U(:,1:T)')./T;   % same as 2:T (because U(:,1) == 0)
    s.UU      = (s.Um1Um1.*T + U(:,T+1)*U(:,T+1)')./(T);   % not T+1, since should be 2:T, only first el = 0.
    s.XU      = (s.Xm1Um1.*T + mu(:,T+1) * U(:,T+1)')./(T);  
end

s.infer = struct('mu', mu); s.infer.P = sigma; s.infer.G = G;   % hack to stop MATLAB making a struct array

% purely for debugging EM
% Q1 = zeros(obj.d.y);
% Q2 = zeros(obj.d.y);
% Xm1_U = zeros(size(s.Xm1_U));
% for tt = 1:T
%     Xm1_U = Xm1_U + mu(:,tt)*U(:,tt+1)';
%     Q1 = Q1 + (mu(:, tt+1) - obj.par.A * mu(:, tt) - obj.par.B * U(:, tt+1))*(mu(:, tt+1) - obj.par.A * mu(:, tt) - obj.par.B * U(:, tt+1))';
%     Q2 = Q2 + (mu(:, tt+1) - obj.par.A * mu(:, tt))*(mu(:, tt+1) - obj.par.A * mu(:, tt))';
%     Q1 = Q1 + sigma{tt+1} - sigma{tt+1}*G{tt}'*obj.par.A' - obj.par.A*G{tt}*sigma{tt+1}' + obj.par.A*sigma{tt}*obj.par.A';
%     Q2 = Q2 + sigma{tt+1} - sigma{tt+1}*G{tt}'*obj.par.A' - obj.par.A*G{tt}*sigma{tt+1}' + obj.par.A*sigma{tt}*obj.par.A';
% end
% s.Q1 = Q1./T;
% s.Q2 = Q2./T;
% s.Xm1_U2 = Xm1_U./T;
end