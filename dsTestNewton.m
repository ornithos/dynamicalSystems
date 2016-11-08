% initial "circular motion" trial
tmp = dynamicalSystem(2, 2, 'x0', 1e-5, 'evolution', [cos(0.03) -sin(0.03); sin(0.03) cos(0.03)], eye(2), ...
        'emission', 0.5*eye(2), 0.2*eye(2), 'data', 50);
tmp = tmp.smoothLinear;
gui.posteriorGaussGUI(tmp);

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
dsNewton  = dynamicalSystem(4, 2, 'x0', 1e-5, 'evolution', transition, Q, 'emission', emission, 5*eye(2), 'data', 100);
dsNewton  = dsNewton.smoothLinear;
dsNewton  = dsNewton.save('original-params');

%%
opts             = struct('maxiter', 1000, 'epsilon', 1e-3);
[dsNewton, llh]  = dsNewton.parameterLearningEM(opts);
dsNewton         = dsNewton.save('learned-all-EM');

dsNewton         = dsNewton.useSavedParameters('original-params');
[dsNewton, llh]  = dsNewton.parameterLearningEM(opts);
dsNewton         = dsNewton.save('learned-all-EM2');

dsNewton         = dsNewton.useSavedParameters('original-params');
[dsNewton, llh]  = dsNewton.parameterLearningEM(opts);
dsNewton         = dsNewton.save('learned-all-EM3');
%%
dsNewtonPlot = dsNewton;
dsNewtonPlot.x = dsNewtonPlot.x(1:2,:);
for ii = 1:dsNewtonPlot.stackptr
    H = dsNewtonPlot.stack{ii,1}.par.H;
    dsNewtonPlot.stack{ii,1}.infer.filter.mu = H * dsNewtonPlot.stack{ii,1}.infer.filter.mu;
    dsNewtonPlot.stack{ii,1}.infer.smooth.mu = H * dsNewtonPlot.stack{ii,1}.infer.smooth.mu;
    for tt = 1:dsNewtonPlot.d.T
        dsNewtonPlot.stack{ii,1}.infer.filter.sigma{tt} = H * dsNewtonPlot.stack{ii,1}.infer.filter.sigma{tt} * H';
        dsNewtonPlot.stack{ii,1}.infer.smooth.sigma{tt} = H * dsNewtonPlot.stack{ii,1}.infer.smooth.sigma{tt} * H';
    end
end
%dsNewtonPlot.filter.mu = H * dsNewtonPlot.filter.mu;
%dsNewtonPlot.smooth.mu = H * dsNewtonPlot.smooth.mu;   % CANNOT OVERRIDE
gui.posteriorGaussGUI(dsNewtonPlot, 'original-params', 'learned-all-EM');