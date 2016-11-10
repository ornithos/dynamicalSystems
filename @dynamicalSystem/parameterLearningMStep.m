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

s = obj.suffStats(struct('bIgnoreHash', true));
% SIGMA = s.SIGMA;
% PHI   = s.PHI;
% B     = s.B;
% D     = s.D;
% C     = s.C;

fOpts = struct('bDoValidation', false, 'bIgnoreHash', true);
% Parameter updates
if isempty(updateOnly)
    if verbose; prevLLH = obj.filter('Kalman',true,[],fOpts); prevLLH  = prevLLH.infer.llh; end
    
    % A and Q updates
    A           = s.C / s.PHI;
    if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''A'' --- Change in LLH: % 5.3f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    Q           = s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A';
    obj.par.A   = A;
    obj.par.Q   = (Q + Q')./2;
    if norm(Q - obj.par.Q)/obj.d.x^2 > 1e-4; warning('Q is not symmetric (%.5f)\n', norm(Q - obj.par.Q)/obj.d.x^2); end
    if verbose; curLLH = obj.filter('Kalman',true,[],fOpts);fprintf('M-Step: ''Q'' --- Change in LLH: % 5.3f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    
    % H and R updates
    H           = s.B / s.SIGMA;
    if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''H'' --- Change in LLH: % 5.3f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    R           = s.D - H*s.B' - s.B*H' + H*s.SIGMA*H';
    obj.par.H   = H;
    obj.par.R   = (R + R')./2;
    if norm(R - obj.par.R)/obj.d.y^2 > 1e-4; warning('R is not symmetric (%.5f)\n', norm(R - obj.par.R)/obj.d.y^2); end
    if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''R'' --- Change in LLH: % 5.3f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end

else
    if verbose; warning('verbose not available when selecting particular M-step updates'); end
    if any(strcmpi(updateOnly, {'A','Q'}))
        A = s.C / s.PHI;
        if any(strcmpi(updateOnly, 'A'))
            obj.par.A = A;
        end
        if any(strcmpi(updateOnly, 'Q'))
            Q           = s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A';
            obj.par.Q   = (Q + Q')./2;
        end
    end
    
    if any(strcmpi(updateOnly, {'H', 'R'}))
        H             = s.B / s.SIGMA;
        if any(strcmpi(updateOnly, 'H'))
            obj.par.H = H;
        end
        if any(strcmpi(updateOnly, 'R'))
            R           = s.D - H*s.B' - s.B*H' + H*s.SIGMA*H';
            obj.par.R   = (R + R')./2;
        end
    end
end


