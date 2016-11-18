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

%% quantites introduced by control inputs
% I'm so sorry for these names...
if any(obj.hasControl)
    u          = [zeros(obj.d.u,1), obj.u];
    u_outer    = u * u';
    u_outer_m  = u_outer - u(:,end)*u(:,end)';
    u_outer_p  = u_outer - u(:,1)*u(:,1)';
end
if obj.hasControl(1)
    Bum     = obj.par.B * u(:,1:end-1) * s.infer.mu(:,1:end-1)'; % subscript is t-1
    Bum     = Bum./(obj.d.T-1);
else
    Bum     = zeros(obj.d.x);
end
if obj.hasControl(2)
    Cum     = obj.par.C * u(:,2:end) * s.infer.mu(:,2:end)'; % subscript is t
    Cum     = Cum./obj.d.T;
else
    Cum     = zeros(obj.d.y,obj.d.x);
end
eps = 5e-4./obj.d.T;
fOpts = struct('bDoValidation', false, 'bIgnoreHash', true);
%% Parameter updates
if isempty(updateOnly)
    if verbose; prevLLH = obj.filter('Kalman',true,[],fOpts); prevLLH  = prevLLH.infer.llh; end
    
    % ============= A and Q updates =======================================
    A           = (s.C - Bum) / s.PHI;
    if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''A'' --- Change in LLH: %5.8f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    Q           = s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A';
    if obj.hasControl(1)
        % --- additions to covariance from control inputs ----.
        Am_m    = obj.par.A * s.infer.mu(:,1:end-1) - s.infer.mu(:,2:end); % sum Am_{t-1} - m_t
        BuAm_m  = obj.par.B * u(:,1:end-1) * Am_m';
        Q       = Q + (obj.par.B*u_outer_m*obj.par.B' + BuAm_m + BuAm_m')./(obj.d.T-1);
    end
    obj.par.A   = A;
    obj.par.Q   = (Q + Q')./2;
    if norm(Q - obj.par.Q)/obj.d.x^2 > 1e-4; warning('Q is not symmetric (%.5f)\n', norm(Q - obj.par.Q)/obj.d.x^2); end
    if verbose; curLLH = obj.filter('Kalman',true,[],fOpts);fprintf('M-Step: ''Q'' --- Change in LLH: %5.8f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    
    % ============= H and R updates =======================================
    H           = (s.B - Cum) / s.SIGMA;
    if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''H'' --- Change in LLH: %5.8f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    R           = s.D - H*s.B' - s.B*H' + H*s.SIGMA*H';
    if obj.hasControl(2)
        % --- additions to covariance from control inputs ----.
        Hm_y    = obj.par.H * s.infer.mu(:,2:end) - obj.y; % sum Hm_{t-1} - y_t
        CuHm_y  = obj.par.C * u(:,2:end) * Hm_y';
        R       = R + (obj.par.C*u_outer_p*obj.par.C' + CuHm_y + CuHm_y')./(obj.d.T);
    end
    obj.par.H   = H;
    obj.par.R   = (R + R')./2;
    if norm(R - obj.par.R)/obj.d.y^2 > 1e-4; warning('R is not symmetric (%.5f)\n', norm(R - obj.par.R)/obj.d.y^2); end
    if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''R'' --- Change in LLH: %5.8f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end

    % ============= B and C updates =======================================
    if obj.hasControl(1)
        B       = (s.infer.mu(:,2:end) - obj.par.A * s.infer.mu(:,1:end-1)) * u(:,1:end-1)';
        B       = B \ (obj.par.uu - obj.u(:,end)*obj.u(:,end)');  % subscript t-1
        obj.par.B = B;
        if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''B'' --- Change in LLH: %5.8f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    end
    if obj.hasControl(2)
        C       = (obj.y - obj.par.H * s.infer.mu(:,2:end)) * u(:,2:end)';
        C       = C \ obj.par.uu;
        obj.par.C = C;
        if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''C'' --- Change in LLH: %5.8f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    end
    
else
    if verbose; prevLLH = obj.filter('Kalman',true,[],fOpts); prevLLH  = prevLLH.infer.llh; end
    % ============= A and Q updates =======================================
    if any(strcmpi(updateOnly, 'A'))
        obj.par.A   = (s.C - Bum) / s.PHI;
        if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''A'' --- Change in LLH: %5.8f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    end
    if any(strcmpi(updateOnly, 'Q'))
        A           = obj.par.A;
        Q           = s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A';
        if obj.hasControl(1)
            % --- additions to covariance from control inputs ----
            Am_m    = A * s.infer.mu(:,1:end-1) - s.infer.mu(:,2:end); % sum Am_{t-1} - m_t
            BuAm_m  = obj.par.B * u(:,1:end-1) * Am_m';
            Q       = Q + (obj.par.B*u_outer_m*obj.par.B' + BuAm_m + BuAm_m')./obj.d.T;
        end

        % -------------------------------------
        % I had a problem with the M-step: here is my alternative Q
        % calc:
%             Q2 = 0;
%             for tt = 2:obj.d.T+1
%                 del = s.infer.mu(:,tt) - A*s.infer.mu(:,tt-1); %- obj.par.B*u(:,tt-1);
%                 Q2  = Q2 + del*del' + s.infer.P{tt} - s.infer.P{tt}*s.infer.G{tt-1}'*A';
%                 Q2  = Q2 - A*s.infer.G{tt-1}*s.infer.P{tt} + A*s.infer.P{tt-1}*A';
%             end
%             Q2 = Q2./(obj.d.T);
%             if norm(Q - Q2) > 1e-3; fprintf('Q diff from Q2! (%.10f)\n', norm(Q - Q2)); keyboard; end
%             Q = Q2;
        % --------------------------------------

%             Q(1:obj.d.x+1:obj.d.x^2) = max(diag(Q), eps);
        obj.par.Q   = (Q + Q')./2 + eps*eye(obj.d.x);

        if verbose; curLLH = obj.filter('Kalman',true,[],fOpts);fprintf('M-Step: ''Q'' --- Change in LLH: %5.8f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    end
    
    % ============= H and R updates =======================================
    if any(strcmpi(updateOnly, 'H'))
        obj.par.H   = (s.B - Cum) / s.SIGMA;
        if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''H'' --- Change in LLH: %5.8f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    end
    if any(strcmpi(updateOnly, 'R'))
        H           = obj.par.H;
        R           = s.D - H*s.B' - s.B*H' + H*s.SIGMA*H';

        % ------ R CAN SOMETIME SHOW A SLIGHT DECREASE (O(1e4)) -------
        % ------ can't find the problem..... --------------------------
        % ------ Also think this is only true when control par B ------
        % ------ is in use....? ---------------------------------------
        % -------------------------------------------------------------

%             if obj.hasControl(2)
%                 % --- additions to covariance from control inputs ----.
%                 Hm_y    = obj.par.H * s.infer.mu(:,2:end) - obj.y; % sum Hm_{t-1} - y_t
%                 CuHm_y  = obj.par.C * u(:,2:end) * Hm_y';
%                 R       = R + (obj.par.C*u_outer_p*obj.par.C' + CuHm_y + CuHm_y')./(obj.d.T);
%             end
%             
%             R2 = 0;
%             for tt = 2:obj.d.T+1
%                 del = (obj.y(:,tt-1) - H*s.infer.mu(:,tt));
%                 if obj.hasControl(2); del = del - obj.par.C*u(:,tt); end
%                 R2  = R2 + del*del' +  H*s.infer.P{tt}*H';
%             end
%             R2 = R2./(obj.d.T);
%             if norm(R - R2) > 1e-9; fprintf('R diff from R2! (%.10f)\n', norm(R - R2)); end
%             Q = R2;
%              fprintf('R diff from R2! (%.15f)\n', norm(R - R2));
%             R(1:obj.d.y+1:obj.d.y^2) = max(diag(R), eps);
        obj.par.R   = (R + R')./2 + eps*eye(obj.d.y);
        if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''R'' --- Change in LLH: %5.8f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    end
    
    % ============= B and C updates =======================================
    if obj.hasControl(1) && any(strcmpi(updateOnly, {'B'}))
        B       = (s.infer.mu(:,2:end) - obj.par.A * s.infer.mu(:,1:end-1)) * u(:,1:end-1)';
        B       = B / (obj.par.uu - obj.u(:,end)*obj.u(:,end)');  % subscript t-1
        obj.par.B = B;
        if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''B'' --- Change in LLH: %5.8f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    end
    if obj.hasControl(2) && any(strcmpi(updateOnly, {'C'}))
        C       = (obj.y - obj.par.H * s.infer.mu(:,2:end)) * u(:,2:end)';
        C       = C / obj.par.uu;
        obj.par.C = C;
        if verbose; curLLH = obj.filter('Kalman',true,[],fOpts); fprintf('M-Step: ''C'' --- Change in LLH: %5.8f\n', curLLH.infer.llh - prevLLH); prevLLH = curLLH.infer.llh; end
    end
end


