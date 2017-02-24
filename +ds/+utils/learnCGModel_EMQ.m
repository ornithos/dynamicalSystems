function [M, exitflag, val] = learnCGModel_EMQ(S1, S2, Vsum, Csum, Q, M_init, M_prev, verbose)
% LEARNCGMODEL Learns the dynamics matrix of a Linear Dynamical System 
% (LDS) from states using the constraint generation algorithm.
%
% (** ADAPTED FROM ORIGINAL FOR:
%     - change metric by Q^{-1}, as per objective function; their EM
%       variant seems to skip out on this.. Perhaps they assume Q = I?
%  ** 2016-12-05  AB ***********************)
%
%  (AB) Note that while it is not really documented, the solver is actually
%  optimising M' not M, and it is converted back in the end.
%
%  Syntax
%  
%    M = learnCGModel_EMQ(S1,S2,Vsum,Csum,Q,M_init,verbose)
%
%  Description
%
%     Implements constraint generation algorithm described in
%     ' A Constraint Generation Approach for Learning Stable Linear
%     Dynamical Systems'
%     Sajid M. Siddiqi, Byron Boots, Geoffrey J.Gordon
%     NIPS 2007
%
%       
%   Given matrices of x_t states and x_{t+1} states, learns a *stable*
%   dynamics matrix that maps x_t to x_{t+1} by starting with the
%   least-squares solution as an unconstrained QP, then repeating:
%   check for stability, use the unstable solution to generate a constraint
%   and add it to the QP. When simulating LB-1 using constraint generation, 
%   constraints are added until the top singular value is 1.

%   Currently, Matlab's 'quadprog' function is used to solve the QP. 
%   CVX (with Sedumi), publicly available
%   optimization software, can also be used. 
%
%
%   M = learnCGModel(S1,S2,simulate_LB1) takes these inputs,
%      S1       - Matrix of state estimates at successive timesteps
%                 x_t,x_{t+1},x_{t+2},...
%      S2       - Matrix of next-state estimates at successive timesteps
%                 x_{t+1},x_{t+2},x_{t+3},...
%      V        - Tensor of state covariance matrices at time t
%      VV       - Tensor of state cross-covariances at successive timesteps
%                
%     
%    and returns
%    M          - dynamics matrix estimate
%
% Authors: Sajid Siddiqi and Byron Boots
%
% Modified:
% 2016    - Konstantinos Georgatzis - for EM variant as detailed in the 
%                                     above paper.
% Nov 16  - Alex Bird               - few minor tweaks for input, verbosity
%                                     and stability.
% Feb 17  - Alex Bird               - Add'l arguments and outputs for
%                                     diagnosing problems, reusing solns

warning off optim:quadprog:SwitchToMedScale;

begintime = clock; 

optimverbose = false;
eps      = 0;
sv_eps   = 0.0004999999999;
dotprint = 5;
maxIter  = 1000;


% we keep track of the following quantities over iterations,
% they can be plotted if desired ...
maxevals = [];
minevals = [];
scores = [];
svals = [];

d = size(M_init,1);

% calculating terms required for the quadratic objective function 
% P,q,r (eqn 4 in paper)

% C = S1*S1';
% C2 = S2*S1';
if isdiag(Q)
    invQ = diag(1./diag(Q));
else
    invQ  = inv(Q);
end

Vsum = Vsum; %assuming always N > dx
tmp = Csum; %assuming always N > dx
% q = tmp(:);
% P = zeros(d^2,d^2);

% for i = 1:d
%     P( (i-1)*d + 1 : (i-1)*d + d,   (i-1)*d + 1 : (i-1)*d + d) = Vsum;
% end
% P = (P +P')/2;    % numerical stability problems..
P    = kron(invQ, Vsum);
P    = (P +P')/2;    % numerical stability problems..
CTQI = Csum*invQ;    % Csum is already C^T
q    = CTQI(:);

%r = trace(S2*S2');
r = 0;

% constraints for QP (initially empty)
G = [];
h = [];

% first M is learned unconstrained
%fprintf('calculating initial M...\n');
M_init = M_init';
M_prev = M_prev';

%M = pinv(S1')*S2';
M = M_init;
lsscore = norm(S1'*M - S2','fro')^2;
if verbose > 1; fprintf('frob score = %.4f\n',lsscore); end

if false
    [tmp_evals,max_e,min_e] = get_eigenthings(M);
    maxevals(end+1) = max_e;
    minevals(end+1) = min_e;

    [u,s,v] = svds(M,1); 
    Morig = M;

    % if we're simulating LB-1, watch for top singular value to
    % dip below 1
    % watch for stability 
    if max_e < 1 && min_e > -1
        if verbose > 1; fprintf('done early! full eigenvals:\n'); end
        M = M';
        return;
    else
        if verbose > 1
            fprintf('initial top eigenval = %.4f\n',max_e);
            fprintf('initial smallest eigenval = %.4f\n',min_e);
        end
    end

scores(end+1) = lsscore;

% calculate constraint based on unstable solution
tmp    = u*v';
ebar   = tmp(:);
G      = [G ; ebar'];
h      = [h ; 1];

svals(end+1) = s;
end

if optimverbose
    options = optimset('maxIter', maxIter);
else
    options = optimset('maxIter', maxIter, 'Display', 'off');
end
% iteratively add constraints, recalculate QP solution

for i = 1:1000
    
    % We want to minimize m'*P*m - 2*q'*m + r, whereas    
    % for quadprog, the objective is to minimize x'*H*x/2  + f'*x
    % with constraints Gx - h <= 0, so have to flip some signs ..
    [m,val,exitflag] = quadprog(2*P,2*(-q),G,h,[],[],[],[],M_init(:), options);

%    CAN ALSO USE CVX INSTEAD OF QUADPROG:
%    val = 0;
%     cvx_begin
%         variable m(d*d)
%         minimize(norm(S1'*reshape(m,d,d)-S2', 'fro'))
%         subject to
%             -G*m + h >= 0;
%     cvx_end
%     

    % adding back the constant term in the objective irrelevant to iqph
    
    score = (val-r);
    scores(end+1) = score;
    
    % reshaping M into the matrix we need
    
    Mprev     = M;
    M         = reshape(m,d,d);  
    
    conscore  = norm(S1'*M - S2','fro')^2;
    
    diffscore = (conscore - lsscore)/lsscore;

    [tmp_evals,max_e,min_e] = get_eigenthings(M);

    maxevals  = [maxevals max_e];
    minevals  = [minevals min_e];

    [u,s,v] = svds(M,1); 

    % if we're simulating LB-1, watch for top singular value to
    % dip below 1, or number of iterations to exceed allowed maximum
    
    % watch for stability 
        
    if( max_e < 1 & min_e > -1 )    % then we're done
        if verbose > 1; fprintf('found M, exiting ...\n'); end
        break;
    else
        if verbose >= 1
            fprintf('.');
            if verbose > 1
                if mod(i,dotprint) == 0
                    fprintf('\n');
                    fprintf('eps: %.3f, diffscore: %.4f, top eval: %.7f, small eval: %.7f\n',eps,diffscore,max_e,min_e);
                end
            end
        end
    end

    scores(end+1) = score;
    svals(end+1) = s;

    % don't need to check largest singular value, because it 
    % must be greater than 1 since the eigenvalue was.
    
    % Add a constraint based on largest singular value
    %        fprintf('top eigval too large, adding constraint ... \n');
    if verbose >=1
        fprintf('.');
        if verbose > 1 && mod(i,dotprint) == 0
            fprintf('\n');
            fprintf('eps: %.3f, diffscore: %.4f, top sval: %.7f\n',eps,diffscore,s);
            fprintf('time so far: %.3f\n', etime(clock,begintime));
        end
    end
    tmp = u*v';
    ebar = tmp(:);
    G = [G ; ebar'];
    h = [h ; 1];
end

maxeig= max(abs(eig(M)));
maxsval = svds(M,1);

if verbose > 1
    fprintf('top eigval after constrained QP: %.7f\n',maxeig);
    fprintf('top sval after constrained QP: %.7f\n',maxsval);
end

% refining the solution: binary search to find boundary of stability region

if true %~simulate_LB1
    
    Mbest = [];
    tol = 0.00001;
    lo = 0;
    hi = 1;

    % interpolating from previous best solution
    Morig = Mprev;
    if true %simulate_LB1 == 0
        while hi-lo > tol
            if verbose >= 1; fprintf(','); end
            alpha = lo + (hi-lo)/2;
            Mbest = (1-alpha)*M + alpha*Morig;
            maxeig = max(abs(eig(Mbest)));
            if (maxeig) > 1
                hi = alpha;
            elseif maxeig < 1
                lo = alpha;
            else    % done!
                break
            end
        end
        Mbest = (1-alpha + tol)*M + (alpha-tol)*Morig;
        M = Mbest;
        maxeig= max(abs(eig(M)));
        maxsval = svds(M,1);

        if verbose > 1
            fprintf('top eigval after binary search: %.7f\n',maxeig);
            fprintf('top sval after binary search: %.7f\n',maxsval);
        end
    end

end

% if initial was actually better (possible due to constraints imposed in
% QP - however it will not be by much)

% returning dynamics matrix in proper orientation

if M_prev(:)'*P*M_prev(:) -2*M_prev(:)'*q < (Mbest(:)'*P*Mbest(:) -2*Mbest(:)'*q) %*1.0001
    M = M_prev';
else
    M = M';
end

end

function [actual_evals,max_e,min_e] = get_eigenthings(M)
% maximum and minimum magnitudes of eigenvalues and corresponding
% eigenvectors, and the actual eigenvalues
[tmp_evecs,tmp_evals] = eig(M);
actual_evals = diag(tmp_evals);
evals = abs(actual_evals);
max_e = max(evals);
min_e = min(evals);

end