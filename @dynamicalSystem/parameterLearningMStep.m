function obj = parameterLearningMStep(obj)
% parameterLearningMStep(obj)
% M-Step of the EM algorithm: maximise likelihood of parameters given the
% current estimate of the latent state x.
% This follows the algorithm in Särkkä.
% The updates are exact and analytic: no iteration is required.


% Check for existence of Smoothed estimates
if ~isfield(obj.posterior, 'smooth') || numel(obj.posterior.smooth.sigma) ~= obj.T
    fprintf('Smoothed estimates not found or incorrect format. Rerunning smoother...\n');
    obj = obj.posteriorFilter(false, false);
end

% concatenate with x0 to provide estimates of x from 0 to T
mu    = [obj.posterior.smooth.x0.mu, obj.posterior.smooth.mu];
sigma = vertcat(obj.posterior.smooth.x0.sigma, obj.posterior.smooth.sigma);
G     = vertcat(obj.posterior.smooth.x0.G, obj.posterior.smooth.G);
y     = obj.y;
T     = obj.T;

% mu    indexes from 0:T     (T+1 elements)
% sigma indexes from 0:T     (T+1 elements)
% G     indexes from 0:(T-1) (T   elements)

A     = obj.A;
H     = obj.H;

%% ALGO
% calculate quantities à la Särkkä
almostAllP = sum(cat(3,sigma{2:T}),3);
almostAllM = mu(:,2:T) * mu(:,2:T)';

SIGMA = (almostAllP + sigma{T+1} + almostAllM + mu(:,T+1) * mu(:,T+1)')./T;
PHI   = (almostAllP + sigma{1} + almostAllM + mu(:,1) * mu(:,1)')./T;
B     = (y * mu(:,2:T+1)')./T;
D     = (y * y')./T;

C     = mu(:,2:T+1) * mu(:,1:T)';
for tt = 1:T
    C = C + sigma{tt+1} * G{tt}';
end
C     = C./T;

% Parameter updates
verbose = false;
if verbose; prevLLH = obj.posteriorFilter(true); end
if verbose; prevLLH  = prevLLH.llh; end
obj.A   = C / PHI;
if verbose; curLLH = obj.posteriorFilter(true); fprintf('M-Step: ''A'' --- Change in LLH: % 5.3f\n', curLLH.llh - prevLLH); prevLLH = curLLH.llh; end
obj.Q   = SIGMA - C*A' - A*C' + A*PHI*A';
obj.Q   = (obj.Q + obj.Q')./2;
if verbose; curLLH = obj.posteriorFilter(true); fprintf('M-Step: ''Q'' --- Change in LLH: % 5.3f\n', curLLH.llh - prevLLH); prevLLH = curLLH.llh; end
obj.H   = B / SIGMA;
if verbose; curLLH = obj.posteriorFilter(true); fprintf('M-Step: ''H'' --- Change in LLH: % 5.3f\n', curLLH.llh - prevLLH); prevLLH = curLLH.llh; end
obj.R   = D - H*B' - B*H' + H*SIGMA*H';
obj.R   = (obj.R + obj.R')./2;
if verbose; curLLH = obj.posteriorFilter(true); fprintf('M-Step: ''R'' --- Change in LLH: % 5.3f\n', curLLH.llh - prevLLH); prevLLH = curLLH.llh; end

end


