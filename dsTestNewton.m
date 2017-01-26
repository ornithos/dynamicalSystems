% initial "circular motion" trial
tmp = ds.dynamicalSystem(2, 2, 'x0', 1e-5, 'evolution', [cos(0.03) -sin(0.03); sin(0.03) cos(0.03)], eye(2), ...
        'emission', 0.5*eye(2), 0.2*eye(2), 'data', 50, struct('verbose', false));
tmp.smooth;

%%
% Newtonian dynamics (4D state space to keep linear system with acc/v)
deltaT              = 0.5;
transition          = eye(4);
transition(1,3)     = deltaT;
transition(2,4)     = deltaT;

emission            = [1 0 0 0; 0 1 0 0];
Q                   = 0.01*eye(4);
Q(3,3)              = 0.2;
Q(4,4)              = 0.2;
% x0                  = struct('mu', zeros(4,1), 'sigma', eye(4)*1e2);
x0                  = [];
dsNewton  = ds.dynamicalSystem(4, 2, 'x0', 1e-5, 'evolution', transition, Q, 'emission', emission, 5*eye(2), 'data', 100);
dsNewton.filter([],true);
dsNewton.smooth;
dsNewton.save('original-params');

%%
opts             = struct('maxiter', 1000, 'epsilon', 1e-3, 'sampleStability', 10, 'multistep', 4);
llh = dsNewton.parameterLearningEM(opts);
dsNewton.filter([],true);
dsNewton.smooth;
dsNewton.save('learned-all-EM1');

% dsNewton.useSavedParameters('original-params');
% dsNewton.parameterLearningEM(opts);
% dsNewton.save('learned-all-EM2');
% 
% dsNewton.useSavedParameters('original-params');
% dsNewton.parameterLearningEM(opts);
% dsNewton.save('learned-all-EM3');
%%

ds.gui.posteriorGaussGUI(dsNewton, 'original-params', 'learned-all-EM1');

%%
% Control inputs test
rng(1)
T      = 100; % must be divisible by 20
assert(mod(T,20)==0, 'T must be divisible by 20');
u_inp  = kron(0.05*rand(2,T/20)-0.025,[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5 2]);
B      = [1, 0; 0, 1; 0.1, 0.1; -0.1, 0.1];
C      = [0,0];
opts   = struct('verbose', true, 'warnings', true);
dsNewton  = ds.dynamicalSystem(4, 2, 'x0', 1e-5, 'evolution', transition, Q, 'emission', emission, 5*eye(2), 'data', 100, 'control', u_inp, B, false, opts);
dsNewton  = dsNewton.smooth;
dsNewton  = dsNewton.save('original-params');

if true
opts      = struct('maxiter', 1000, 'epsilon', 1e-3, 'dbg', true);
dsNewton  = dsNewton.parameterLearningEM(opts);
dsNewton  = dsNewton.save('learned-all-EM');
ds.gui.posteriorGaussGUI(dsNewton, 'initialised', 'learned-all-EM');
end

%%
% Online learning test
deltaT              = 0.5;
transition          = eye(4);
transition(1,3)     = deltaT;
transition(2,4)     = deltaT;

emission            = [1 0 0 0; 0 1 0 0];
Q                   = 0.01*eye(4);
Q(3,3)              = 0.2;
Q(4,4)              = 0.2;
% x0                  = struct('mu', zeros(4,1), 'sigma', eye(4)*1e2);
x0                  = [];
dsNewton  = ds.dynamicalSystem(4, 2, 'x0', 1e-5, 'evolution', transition, Q, 'emission', emission, 5*eye(2), 'data', 100);
dsNewton  = dsNewton.smooth;
dsNewton  = dsNewton.save('original-params');
dsNewton.par.A = eye(size(tran

