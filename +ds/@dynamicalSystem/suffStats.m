function s = suffStats(obj, opts)
% s = suffStats(obj, opts)
% Calculates the various quantities used in Särkkä's parameter estimation
% equations. They are typically pairwise expectations.
% opts: verbose, bIgnoreHash, bDoValidation
%
% OUTPUTS:
%  - SIGMA = <x_t, x_t>
%  - PHI   = <x_{t-1}, x_{t-1}>
%  - C     = <x_t, x_{t-1}> 
%  - B     = <y_t, x_t>          only for t in tau
%  - D     = <y_t, y_t>          only for t in tau
% (-SIGMEV = <x_t, x_t>          only for t in tau  <- needed for R in presence of miss. vals)
%
%
%  - XU    = <x_t, u_t>
%  - Xm1_U = <x_{t-1}, u_t>
%  - UU    = <u_t, u_t>

if nargin < 2 || isempty(opts)
    opts = struct;
end

optsDefault  = struct('verbose', true, 'bIgnoreHash', false, 'anneal', 1);
opts         = utils.base.parse_argumentlist(opts, optsDefault);
    
% Check for existence of Smoothed estimates
if ~opts.bIgnoreHash && obj.parametersChanged
    if opts.verbose; fprintf('Filter not run for current params. Rerunning filter/smoother...\n'); end
    if isempty(obj.infer.sType); fType = 'Linear'; else, fType = obj.infer.sType; end
    obj.filter(fType, false, [], opts);
    obj.smooth(fType, [], opts);
end
if isempty(obj.infer.sType)
    if opts.verbose; fprintf('Smoother not run for current params. Running smoother...\n'); end
    obj.smooth('Linear', [], opts);
end

% concatenate with x0 to provide estimates of x from 0 to T
mu    = [obj.infer.smooth.x0.mu, obj.infer.smooth.mu];
sigma = vertcat(obj.infer.smooth.x0.sigma, obj.infer.smooth.sigma);
G     = vertcat(obj.infer.smooth.x0.G, obj.infer.smooth.G);

% deal with missing values in y.
y     = obj.impute_y('smooth', true, 'bIgnoreHash', true);
yActv = ~all(isnan(obj.y),1);
T     = find(yActv, 1, 'last');
yActv = yActv(1:T);
y     = y(:,1:T);
y(:,~yActv) = 0;   % most suff stats are multiplicative, so it is convenient to zero all NaN values.


% deterministic annealing of posterior:
if opts.anneal < 1
    cellfun(@(x) x.*opts.anneal, sigma, 'UniformOutput',0);
    cellfun(@(x) x.*opts.anneal, G, 'UniformOutput',0);
end

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

% transition statistics
s.SIGMA = (almostAllP + sigma{T+1} + almostAllM + mu(:,T+1) * mu(:,T+1)')./T;
s.PHI   = (almostAllP + sigma{1} + almostAllM + mu(:,1) * mu(:,1)')./T;

s.C     = mu(:,2:T+1) * mu(:,1:T)';
for tt = 1:T
    s.C = s.C + sigma{tt+1} * G{tt}';
end
s.C     = s.C./T;

% evolution statistics
% yActv   = ~all(isnan(obj.y),1);   ---> see above
ySemi   = any(isnan(obj.y(:,1:T)), 1);
Ty      = sum(yActv);
s.B     = (y * mu(:,2:T+1)')./Ty;   % only using non all-NaNs: since 0 for all NaNs (see above), and T is last non NaN.
s.D     = (y * y')./Ty;           % only using non all-NaNs


% ------------- MISSING VALUES -------------------------------
% remove x_t terms which are not active
s.SIGMEV= s.SIGMA.*T;
for tt = find(~yActv)
    s.SIGMEV = s.SIGMEV - sigma{tt+1} - mu(:,tt+1) * mu(:,tt+1)';   % + 1 because of x0
end
s.SIGMEV = s.SIGMEV./Ty;
% add covariance into y where relevant
for tt = find(yActv & ySemi)
    mask = isnan(obj.y(:,tt));
    Hu   = obj.par.H(mask,:);
    Ho   = obj.par.H(~mask,:);
    Ru   = obj.par.R(mask, mask);
    Ro   = obj.par.R(~mask, ~mask);
    % D = < y_t, y_t >
    adj  = (Hu * sigma{tt+1} * Hu' + Ru)./Ty;
    s.D(mask,mask) = s.D(mask,mask) + adj;
    % B = < y_t, x_t >
    % THIS IS TOTALLY WRONG. NEED THE JOINT DISTRIBUTION OF y_t, x_t | y_{1:T}
    % WILL GIVE US MEAN
%     y_t        = y(:,tt);
%     y_t(~mask) = 0;   % zero out the give (delta fn) elements of y as do not want to remove these.
%     s.B        = s.B - (y_t* mu(:,tt+1)')./Ty;  
%     HoP        = Ho * sigma{tt+1};
%     adj        = ((Hu * HoP') / ( Ho * HoP' + Ro)) * (HoP * Hu');
    adj        = (Hu * sigma{tt})./Ty;
    s.B(mask, :) = s.B(mask, :) + adj;
end


% ----- Control quantities -----------------------------------
if obj.hasControl(1)
    s.X_Um1   = (mu(:,2:T+1) * U(:,1:T)')./T;
    s.Xm1_U   = (mu(:,1:T)   * U(:,2:T+1)')./T;
    s.U_Um1   = (U(:,2:T+1)  * U(:,1:T)')./T;
end

if obj.hasControl(2)
    s.YU      = (y(:,yActv) * U(:, find(yActv)+1)')./T;
end

if any(obj.hasControl)
    s.Um1Um1  = (U(:,1:T) * U(:,1:T)')./T;
    s.Xm1Um1  = (mu(:,1:T) * U(:,1:T)')./T;   % same as 2:T (because U(:,1) == 0)
    s.UU      = (s.Um1Um1.*T + U(:,T+1)*U(:,T+1)')./(T);   % not T+1, since should be 2:T, only first el = 0.
    s.XU      = (s.Xm1Um1.*T + mu(:,T+1) * U(:,T+1)')./(T);  
end

s.infer = struct('mu', mu); s.infer.P = sigma; s.infer.G = G;   % hack to stop MATLAB making a struct array

% additional single variable means
s.Ymu    = sum(y, 2)./T;
s.Umu    = sum(U(: , 2:T+1), 2)./T;
s.Xmu    = sum(mu(:, 2:T+1), 2)./T;
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