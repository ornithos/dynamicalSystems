function obj = parameterLearningMStep(obj, verbose, updateOnly)
% parameterLearningMStep(obj)
% M-Step of the EM algorithm: maximise likelihood of parameters given the
% current estimate of the latent state x.
% This follows the algorithm in Särkkä.
% The updates are exact and analytic: no iteration is required.

if nargin < 3
    updateOnly = [];
end

if nargin < 2 || isempty(verbose)
    verbose = false;
end

% Check for existence of Smoothed estimates
if obj.infer.fpHash ~= obj.parameterHash
    fprintf('Filter not run or parameters changed. Rerunning filter...\n');
    obj = obj.filterKalman(false, false);
    obj = obj.smoothLinear(false, false);
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

A     = obj.par.A;
H     = obj.par.H;

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
if isempty(updateOnly)
    if verbose; prevLLH = obj.filterKalman(true); prevLLH  = prevLLH.infer.llh; end
    
    % A and Q updates
    A           = C / PHI;
    if verbose; curLLH = obj.filterKalman(true); fprintf('M-Step: ''A'' --- Change in LLH: % 5.3f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    Q           = SIGMA - C*A' - A*C' + A*PHI*A';
    obj.par.A   = A;
    obj.par.Q   = (Q + Q')./2;
    if norm(Q - obj.par.Q)/obj.d.x^2 > 1e-4; warning('Q is not symmetric (%.5f)\n', norm(Q - obj.par.Q)/obj.d.x^2); end
    if verbose; curLLH = obj.filterKalman(true); fprintf('M-Step: ''Q'' --- Change in LLH: % 5.3f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    
    % H and R updates
    H           = B / SIGMA;
    if verbose; curLLH = obj.filterKalman(true); fprintf('M-Step: ''H'' --- Change in LLH: % 5.3f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    R           = D - H*B' - B*H' + H*SIGMA*H';
    obj.par.H   = H;
    obj.par.R   = (R + R')./2;
    if norm(R - obj.par.R)/obj.d.y^2 > 1e-4; warning('R is not symmetric (%.5f)\n', norm(R - obj.par.R)/obj.d.y^2); end
    if verbose; curLLH = obj.filterKalman(true); fprintf('M-Step: ''R'' --- Change in LLH: % 5.3f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end

else
    if verbose; warning('verbose not available when selecting particular M-step updates'); end
    if any(strcmpi(updateOnly, {'A','Q'}))
        A = C / PHI;
        if any(strcmpi(updateOnly, 'A'))
            obj.par.A = A;
        end
        if any(strcmpi(updateOnly, 'Q'))
            Q           = SIGMA - C*A' - A*C' + A*PHI*A';
            obj.par.Q   = (Q + Q')./2;
        end
    end
    
    if any(strcmpi(updateOnly, {'H', 'R'}))
        H             = B / SIGMA;
        if any(strcmpi(updateOnly, 'H'))
            obj.par.H = H;
        end
        if any(strcmpi(updateOnly, 'R'))
            R           = D - H*B' - B*H' + H*SIGMA*H';
            obj.par.R   = (R + R')./2;
        end
    end
end


